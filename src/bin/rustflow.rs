//! Native JSON steering for graph evaluation and progressive physical transport.
use serde::Deserialize;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    sync::Arc,
};
use symbolica::prelude::*;
use symbolica_amflow::{
    algebraic::{AlgebraicKinematicSystem, CanonicalAlgebraicSystem, SquareRoot},
    kinematics::KinematicSystem,
    physical_transport::{BoundaryAttemptOutcome, PhysicalResult},
    transport_cache::*,
    *,
};

type CliResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

#[derive(Deserialize)]
#[serde(default, deny_unknown_fields)]
struct Settings {
    digits: u32,
    guard_digits: u32,
    series_order: usize,
    max_steps: usize,
    max_precision_attempts: usize,
    max_boundary_attempts: usize,
    workers: usize,
    dimension: i64,
    recursion: String,
    prescription: String,
    refine_basis: bool,
    skip_reduction: bool,
    sampled_reduction: bool,
    symbolic_cache: Option<PathBuf>,
    progress: bool,
}
impl Default for Settings {
    fn default() -> Self {
        let f = FlowOptions::default();
        Self {
            digits: f.digits,
            guard_digits: f.guard_digits,
            series_order: f.series_order,
            max_steps: f.max_steps,
            max_precision_attempts: f.max_precision_attempts,
            max_boundary_attempts: f.max_boundary_attempts,
            workers: f.workers,
            dimension: f.dimension,
            recursion: "amf".into(),
            prescription: "+i0".into(),
            refine_basis: f.refine_basis,
            skip_reduction: f.skip_reduction,
            sampled_reduction: f.sampled_reduction,
            symbolic_cache: None,
            progress: false,
        }
    }
}
impl Settings {
    fn options(&self, directory: &Path) -> CliResult<FlowOptions> {
        let options = FlowOptions {
            digits: self.digits,
            guard_digits: self.guard_digits,
            series_order: self.series_order,
            max_steps: self.max_steps,
            max_precision_attempts: self.max_precision_attempts,
            max_boundary_attempts: self.max_boundary_attempts,
            workers: self.workers,
            dimension: self.dimension,
            recursion: match self.recursion.as_str() {
                "amf" => RecursionMode::Amf,
                "ft" => RecursionMode::Ft,
                _ => return Err("recursion must be amf or ft".into()),
            },
            prescription: match self.prescription.as_str() {
                "+i0" => Prescription::PlusI0,
                "-i0" => Prescription::MinusI0,
                _ => return Err("prescription must be +i0 or -i0".into()),
            },
            refine_basis: self.refine_basis,
            skip_reduction: self.skip_reduction,
            sampled_reduction: self.sampled_reduction,
            cache_directory: self.symbolic_cache.as_ref().map(|p| directory.join(p)),
            ..Default::default()
        };
        options.validate()?;
        Ok(options)
    }
    fn context(&self) -> RunContext {
        RunContext {
            progress: self.progress.then(|| {
                Arc::new(|event| eprintln!("{event:?}")) as Arc<dyn Fn(Progress) + Send + Sync>
            }),
            ..Default::default()
        }
    }
}

fn namespace() -> String {
    "feynkit_graph".into()
}
fn epsilon() -> String {
    "eps".into()
}
fn tensor_dimension() -> String {
    "D".into()
}
fn one() -> String {
    "1".into()
}
fn partial_fraction_limit() -> usize {
    10000
}
fn parse(text: &str, namespace: &str) -> CliResult<Atom> {
    Atom::parse(text, namespace.to_owned(), Default::default())
        .map_err(|e| format!("expression {text:?}: {e}").into())
}
fn symbol(text: &str, namespace: &str) -> CliResult<Symbol> {
    match parse(text, namespace)?.as_view() {
        AtomView::Var(v) => Ok(v.get_symbol()),
        _ => Err(format!("expected a symbol, got {text:?}").into()),
    }
}
fn coordinates(
    values: &BTreeMap<String, String>,
    namespace: &str,
) -> CliResult<BTreeMap<Symbol, Atom>> {
    values
        .iter()
        .map(|(k, v)| Ok((symbol(k, namespace)?, parse(v, namespace)?)))
        .collect()
}
fn complex_json(value: &ComplexFloat) -> Value {
    json!({"real":value.re.as_raw().to_string(),"imaginary":value.im.as_raw().to_string()})
}
fn point_json(point: &CachedPoint) -> CliResult<Value> {
    Ok(json!(
        point
            .restart_coordinates()?
            .iter()
            .map(|(k, v)| (Atom::var(*k).to_canonical_string(), v.to_canonical_string()))
            .collect::<BTreeMap<_, _>>()
    ))
}
fn version(schema: u32) -> CliResult<()> {
    if schema != 1 {
        return Err(format!("unsupported steering schema {schema}; expected 1").into());
    }
    Ok(())
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ScalarProduct {
    left: String,
    right: String,
    value: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct GraphRequest {
    schema_version: u32,
    model: PathBuf,
    diagram: PathBuf,
    #[serde(default = "namespace")]
    namespace: String,
    #[serde(default = "epsilon")]
    epsilon: String,
    #[serde(default = "tensor_dimension")]
    tensor_dimension: String,
    #[serde(default)]
    scalar_products: Vec<ScalarProduct>,
    #[serde(default)]
    substitutions: BTreeMap<String, String>,
    #[serde(default)]
    edge_powers: BTreeMap<usize, i16>,
    #[serde(default)]
    last_epsilon_power: i32,
    #[serde(default = "partial_fraction_limit")]
    max_partial_fraction_states: usize,
    #[serde(default)]
    options: Settings,
}

fn graph(request: GraphRequest, directory: &Path) -> CliResult<Value> {
    use feynkit_graph::EdgeId;
    use feynkit_kinematics::Kinematics;
    use feynkit_model::Model;
    use symbolica_amflow::hepkit::GraphIntegral;
    version(request.schema_version)?;
    let options = request.options.options(directory)?;
    let context = request.options.context();
    let ns = &request.namespace;
    // Native tensor structures require a symbolic dimension until contraction.
    let mut kinematics = Kinematics::in_dimension(&parse(&request.tensor_dimension, ns)?)?;
    for product in request.scalar_products {
        kinematics = kinematics.with_scalar_product(
            &parse(&product.left, ns)?,
            &parse(&product.right, ns)?,
            parse(&product.value, ns)?,
        )?;
    }
    let point = KinematicPoint(
        request
            .substitutions
            .iter()
            .map(|(k, v)| Ok((parse(k, ns)?, parse(v, ns)?)))
            .collect::<CliResult<_>>()?,
    );
    let model = Arc::new(Model::from_json(&fs::read_to_string(
        directory.join(request.model),
    )?)?);
    let diagram = fs::read_to_string(directory.join(request.diagram))?;
    let input = GraphIntegral::from_dot(model, &diagram, &kinematics)?.with_powers(
        &request
            .edge_powers
            .iter()
            .map(|(edge, power)| (EdgeId(*edge), *power))
            .collect(),
    )?;
    let groups = input.integral_groups(
        &point,
        symbol(&request.epsilon, ns)?,
        options.dimension,
        request.max_partial_fraction_states,
        &context,
    )?;
    let result = solve_integral_combinations(
        &groups,
        &KinematicPoint::default(),
        request.last_epsilon_power,
        &options,
        &RustRedBackend::default(),
        &context,
    )?;
    Ok(json!({
        "schema_version":1,"operation":"graph","verified_digits":result.verified_digits,
        "working_bits":result.working_bits,"samples":result.samples,
        "validation_samples":result.validation_samples,"refinements":result.refinements,
        "coefficients":result.coefficients.iter().map(|(k,v)| (k.to_string(),complex_json(v))).collect::<BTreeMap<_,_>>(),
        "comparison_errors":result.comparison_errors.iter().map(|(k,v)| (k.to_string(),v.as_raw().to_string())).collect::<BTreeMap<_,_>>()
    }))
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DecimalComplex {
    real: String,
    imaginary: String,
}
#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum SheetChoice {
    Principal,
    Opposite,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct DetailedDestination {
    coordinates: BTreeMap<String, String>,
    #[serde(default)]
    root_germ: Option<BTreeMap<String, SheetChoice>>,
}

#[derive(Deserialize)]
#[serde(untagged)]
enum Destination {
    Detailed(DetailedDestination),
    // Preserve the original rational transport coordinate-map format.
    Coordinates(BTreeMap<String, String>),
}
impl Destination {
    fn parts(
        &self,
    ) -> (
        &BTreeMap<String, String>,
        Option<&BTreeMap<String, SheetChoice>>,
    ) {
        match self {
            Self::Detailed(value) => (&value.coordinates, value.root_germ.as_ref()),
            Self::Coordinates(value) => (value, None),
        }
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Seed {
    coordinates: BTreeMap<String, String>,
    #[serde(default)]
    root_germ: Option<BTreeMap<String, SheetChoice>>,
    working_digits: u32,
    verified_digits: u32,
    coefficients: Vec<Vec<DecimalComplex>>,
    absolute_errors: Vec<Vec<String>>,
    provenance: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct CanonicalInput {
    variables: Vec<String>,
    letters: Vec<String>,
    matrices: Vec<Vec<Vec<String>>>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TransportRequest {
    schema_version: u32,
    #[serde(default = "namespace")]
    namespace: String,
    #[serde(default = "epsilon")]
    epsilon: String,
    #[serde(default)]
    derivatives: Option<BTreeMap<String, Vec<Vec<String>>>>,
    #[serde(default)]
    canonical: Option<CanonicalInput>,
    #[serde(default)]
    roots: BTreeMap<String, String>,
    basis: Vec<String>,
    #[serde(default = "one")]
    normalization: String,
    branch_domain: String,
    #[serde(default)]
    nonzero_conditions: Vec<String>,
    leading_epsilon_power: i32,
    last_epsilon_power: i32,
    cache_directory: PathBuf,
    #[serde(default)]
    seeds: Vec<Seed>,
    destinations: Vec<Destination>,
    #[serde(default)]
    distance_scales: BTreeMap<String, String>,
    #[serde(default)]
    options: Settings,
}
fn transport(request: TransportRequest, directory: &Path) -> CliResult<Value> {
    version(request.schema_version)?;
    let options = request.options.options(directory)?;
    let context = request.options.context();
    let ns = &request.namespace;
    if request.derivatives.is_some() == request.canonical.is_some() {
        return Err("transport requires exactly one of derivatives or canonical".into());
    }
    let roots = request
        .roots
        .iter()
        .map(|(name, radicand)| {
            Ok(SquareRoot {
                symbol: symbol(name, ns)?,
                radicand: parse(radicand, ns)?,
            })
        })
        .collect::<CliResult<Vec<_>>>()?;
    let basis = request
        .basis
        .iter()
        .map(|a| parse(a, ns))
        .collect::<CliResult<Vec<_>>>()?;
    let normalization = parse(&request.normalization, ns)?;
    let conditions = request
        .nonzero_conditions
        .iter()
        .map(|a| parse(a, ns))
        .collect::<CliResult<Vec<_>>>()?;
    let range = EpsilonRange::new(request.leading_epsilon_power, request.last_epsilon_power)?;
    if let Some(canonical) = &request.canonical {
        let system = CanonicalAlgebraicSystem::new(
            symbol(&request.epsilon, ns)?,
            &canonical
                .variables
                .iter()
                .map(|s| symbol(s, ns))
                .collect::<CliResult<Vec<_>>>()?,
            &canonical
                .letters
                .iter()
                .map(|a| parse(a, ns))
                .collect::<CliResult<Vec<_>>>()?,
            &canonical
                .matrices
                .iter()
                .map(|matrix| {
                    matrix
                        .iter()
                        .map(|row| row.iter().map(|a| parse(a, ns)).collect())
                        .collect()
                })
                .collect::<CliResult<Vec<_>>>()?,
            roots,
        )?;
        let flow = RustFlow::with_canonical_conditions(
            system,
            &basis,
            &normalization,
            options.prescription,
            &request.branch_domain,
            &conditions,
        )?;
        let mut result = execute_transport(
            &request,
            directory,
            flow.identity(),
            range,
            |cache, point, policy| {
                Ok(flow.evaluate_to(
                    cache,
                    &point.restart_coordinates()?,
                    point.root_germ(),
                    range,
                    &options,
                    &context,
                    policy,
                )?)
            },
        )?;
        result["representation"] = "canonical".into();
        result["identity"] = flow.identity().key().into();
        return Ok(result);
    }
    let derivatives = request.derivatives.as_ref().ok_or("missing derivatives")?;
    let system = KinematicSystem {
        epsilon: symbol(&request.epsilon, ns)?,
        derivatives: derivatives
            .iter()
            .map(|(name, matrix)| {
                Ok((
                    symbol(name, ns)?,
                    matrix
                        .iter()
                        .map(|row| row.iter().map(|a| parse(a, ns)).collect())
                        .collect::<CliResult<_>>()?,
                ))
            })
            .collect::<CliResult<_>>()?,
    };
    if request.roots.is_empty() {
        let flow = RustFlow::with_conditions(
            system,
            &basis,
            &normalization,
            options.prescription,
            &request.branch_domain,
            &conditions,
        )?;
        execute_transport(
            &request,
            directory,
            flow.identity(),
            range,
            |cache, point, policy| {
                Ok(flow.evaluate_to(
                    cache,
                    &point.restart_coordinates()?,
                    range,
                    &options,
                    &context,
                    policy,
                )?)
            },
        )
    } else {
        let flow = RustFlow::with_algebraic_conditions(
            AlgebraicKinematicSystem {
                epsilon: system.epsilon,
                derivatives: system.derivatives,
                roots,
            },
            &basis,
            &normalization,
            options.prescription,
            &request.branch_domain,
            &conditions,
        )?;
        execute_transport(
            &request,
            directory,
            flow.identity(),
            range,
            |cache, point, policy| {
                Ok(flow.evaluate_to(
                    cache,
                    &point.restart_coordinates()?,
                    point.root_germ().ok_or("missing destination root_germ")?,
                    range,
                    &options,
                    &context,
                    policy,
                )?)
            },
        )
    }
}

fn transport_point(
    values: &BTreeMap<String, String>,
    choices: Option<&BTreeMap<String, SheetChoice>>,
    identity: &BoundaryIdentity,
    namespace: &str,
) -> CliResult<CachedPoint> {
    let point = CachedPoint::Exact(coordinates(values, namespace)?);
    if identity.roots().is_empty() {
        if choices.is_some() {
            return Err("root_germ requires a nonempty roots registry".into());
        }
        return Ok(point);
    }
    let choices =
        choices.ok_or("each algebraic seed and destination requires an explicit root_germ")?;
    let mut sheets = BTreeMap::new();
    for (name, choice) in choices {
        let root = symbol(name, namespace)?;
        let sheet = match choice {
            SheetChoice::Principal => RootSheet::Principal,
            SheetChoice::Opposite => RootSheet::Opposite,
        };
        if sheets.insert(root, sheet).is_some() {
            return Err(format!("root_germ names the same root twice: {name}").into());
        }
    }
    if sheets.len() != identity.roots().len()
        || identity
            .roots()
            .iter()
            .any(|r| !sheets.contains_key(&r.symbol))
    {
        return Err(
            "root_germ must name every registered root exactly once and no other symbol".into(),
        );
    }
    Ok(point.with_root_germ(RootGerm { sheets })?)
}

fn germ_json(point: &CachedPoint) -> Option<Value> {
    point.root_germ().map(|germ| {
        json!(
            germ.sheets
                .iter()
                .map(|(&s, sheet)| (
                    Atom::var(s).to_canonical_string(),
                    match sheet {
                        RootSheet::Principal => "principal",
                        RootSheet::Opposite => "opposite",
                    },
                ))
                .collect::<BTreeMap<_, _>>()
        )
    })
}

// Both native connection kinds share parsing, bank I/O, progressive insertion,
// policy and result serialization. Numerical transport stays in RustFlow.
fn execute_transport(
    request: &TransportRequest,
    directory: &Path,
    identity: &BoundaryIdentity,
    range: EpsilonRange,
    mut evaluate: impl FnMut(
        &mut RustFlowCache,
        &CachedPoint,
        &dyn TransportCost,
    ) -> CliResult<PhysicalResult>,
) -> CliResult<Value> {
    let ns = &request.namespace;
    let cache_directory = directory.join(&request.cache_directory);

    let mut cache = RustFlowCache::load(&cache_directory)?;
    let loaded = cache.len();
    let mut seeds = Vec::new();
    for seed in &request.seeds {
        let p = Precision::decimal(seed.working_digits)?;
        seeds.push(CachedBoundary {
            identity: identity.clone(),
            point: transport_point(&seed.coordinates, seed.root_germ.as_ref(), identity, ns)?,
            kind: PointKind::Physical,
            range,
            coefficients: seed
                .coefficients
                .iter()
                .map(|row| row.iter().map(|c| p.parse(&c.real, &c.imaginary)).collect())
                .collect::<symbolica_amflow::Result<_>>()?,
            accuracy: BoundaryAccuracy::supplied(
                seed.verified_digits,
                p.bits,
                seed.absolute_errors
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|a| {
                                Float::parse(a, Some(p.bits))
                                    .map_err(|e| Error::InvalidInput(format!("seed error: {e}")))
                            })
                            .collect()
                    })
                    .collect::<symbolica_amflow::Result<_>>()?,
                &seed.provenance,
            )?,
        });
    }
    cache.insert_many(seeds)?;
    let policy = ScaledDistance {
        scales: coordinates(&request.distance_scales, ns)?,
        // The request supplies a common branch domain. The physical orchestrator
        // also rejects singularities on each selected straight parameter path.
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let mut results = Vec::new();
    for destination in &request.destinations {
        let (coordinates, germ) = destination.parts();
        let point = transport_point(coordinates, germ, identity, ns)?;
        let result = evaluate(&mut cache, &point, &policy)?;
        cache.save(&cache_directory)?;
        let mut output = json!({
            "coordinates":point_json(&result.boundary.point)?,
            "starting_point":point_json(&result.starting_point)?,
            "verified_digits":result.boundary.accuracy.verified_digits(),
            "working_bits":result.boundary.accuracy.working_bits(),
            "steps":result.transport.as_ref().map_or(0,|t| t.diagnostics.steps),
            "inserted_points":result.inserted_points,
            "coefficients":result.boundary.coefficients.iter().map(|row| row.iter().map(complex_json).collect::<Vec<_>>()).collect::<Vec<_>>(),
            "absolute_errors":result.boundary.accuracy.comparison_errors().iter().map(|row| row.iter().map(|a| a.as_raw().to_string()).collect::<Vec<_>>()).collect::<Vec<_>>()
        });
        if result.boundary_attempts.len() > 1 {
            output["boundary_attempts"] = Value::Array(result.boundary_attempts.iter().map(|attempt| {
                let mut entry = json!({
                    "starting_point":point_json(&attempt.starting_point)?,
                    "cost":attempt.cost.as_raw().to_string(),
                    "source_verified_digits":attempt.source_verified_digits,
                    "outcome":match &attempt.outcome {
                        BoundaryAttemptOutcome::Accepted => json!({"status":"accepted"}),
                        BoundaryAttemptOutcome::AccuracyRejected { message } => json!({"status":"accuracy_rejected","message":message}),
                    },
                });
                if let Some(germ) = germ_json(&attempt.starting_point) {entry["root_germ"] = germ;}
                Ok(entry)
            }).collect::<CliResult<Vec<_>>>()?);
        }
        if let Some(germ) = germ_json(&result.boundary.point) {
            output["root_germ"] = germ;
            output["starting_root_germ"] =
                germ_json(&result.starting_point).ok_or("missing selected source root_germ")?;
        }
        results.push(output);
    }
    cache.save(&cache_directory)?;
    Ok(
        json!({"schema_version":1,"operation":"transport","loaded_boundaries":loaded,"retained_boundaries":cache.len(),"results":results}),
    )
}

fn run() -> CliResult<()> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if arguments.is_empty() || matches!(arguments[0].as_str(), "--help" | "-h" | "help") {
        println!(
            "RustFlow\n\nUsage: rustflow graph INPUT.json\n       rustflow transport INPUT.json\n\nGraph input references native HEPKit model JSON and DOT files.\nTransport input supplies differential equations and a binary cache directory.\nPaths are relative to INPUT.json; results are JSON on stdout.\nSee docs/cli.md for exact-input and accuracy conventions."
        );
        return Ok(());
    }
    if arguments == ["--version"] {
        println!("RustFlow {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if arguments.len() != 2 {
        return Err("expected a command and one steering JSON file; use --help".into());
    }
    let path = Path::new(&arguments[1]);
    let directory = path.parent().unwrap_or_else(|| Path::new("."));
    let text = fs::read_to_string(path)?;
    let output = match arguments[0].as_str() {
        "graph" => graph(serde_json::from_str(&text)?, directory)?,
        "transport" => transport(serde_json::from_str(&text)?, directory)?,
        _ => return Err("unknown command; expected graph or transport".into()),
    };
    println!("{}", serde_json::to_string_pretty(&output)?);
    Ok(())
}
fn main() {
    if let Err(error) = run() {
        eprintln!("RustFlow: {error}");
        std::process::exit(1);
    }
}
