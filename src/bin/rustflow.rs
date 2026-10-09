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
    contour::PolynomialPrescription,
    kinematics::KinematicSystem,
    physical_transport::{BoundaryAttemptOutcome, PhysicalResult, PhysicalRoute},
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
    pade_degree: Option<usize>,
    residual_arithmetic: String,
    boundary_error_strategy: String,
    max_steps: usize,
    step_size_strategy: String,
    local_coordinate: String,
    max_precision_attempts: usize,
    max_boundary_attempts: usize,
    workers: usize,
    dimension: i64,
    recursion: String,
    prescription: String,
    mass_mode: Option<String>,
    deformed_propagator_slots: Option<Vec<usize>>,
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
            pade_degree: None,
            residual_arithmetic: "ball".into(),
            boundary_error_strategy: "automatic".into(),
            max_steps: f.max_steps,
            step_size_strategy: "halving".into(),
            local_coordinate: "identity".into(),
            max_precision_attempts: f.max_precision_attempts,
            max_boundary_attempts: f.max_boundary_attempts,
            workers: f.workers,
            dimension: f.dimension,
            recursion: "amf".into(),
            prescription: "+i0".into(),
            mass_mode: None,
            deformed_propagator_slots: None,
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
            pade: self.pade_degree.map(|degree| PadeOptions {
                degree,
                ..Default::default()
            }),
            residual_arithmetic: self.residual_arithmetic.parse()?,
            boundary_error_strategy: self.boundary_error_strategy.parse()?,
            max_steps: self.max_steps,
            step_size_strategy: match self.step_size_strategy.as_str() {
                "halving" => StepSizeStrategy::Halving,
                "bracketed" => StepSizeStrategy::Bracketed,
                _ => return Err("step_size_strategy must be halving or bracketed".into()),
            },
            local_coordinate: match self.local_coordinate.as_str() {
                "identity" => LocalCoordinate::Identity,
                "balanced_mobius" => LocalCoordinate::BalancedMobius,
                _ => return Err("local_coordinate must be identity or balanced_mobius".into()),
            },
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
            mass_mode: MassMode::from_selection(
                self.mass_mode.as_deref(),
                self.deformed_propagator_slots.clone(),
            )?,
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
struct CutSelection {
    index: usize,
    future_channel: Vec<String>,
    loop_prescriptions: Vec<String>,
    /// Present for finite samples; omitted for independently refined Laurent fitting.
    #[serde(default)]
    epsilon_samples: Option<Vec<String>>,
}

fn exact_rational(value: &str, namespace: &str) -> CliResult<Rational> {
    if let AtomView::Num(n) = parse(value, namespace)?.as_view()
        && let symbolica::coefficient::Coefficient::Complex(c) = n.get_coeff_view().to_owned()
        && c.im.is_zero()
    {
        return Ok(c.re);
    }
    Err(format!("expected an exact real rational, got {value:?}").into())
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
    cut: Option<CutSelection>,
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
    let backend = RustRedBackend::default();
    let (result, operation, conditions) = if let Some(cut) = request.cut {
        let channel = cuts::FutureTimelikeChannel {
            external: cut
                .future_channel
                .iter()
                .map(|value| exact_rational(value, ns))
                .collect::<CliResult<_>>()?,
        };
        let prescriptions = cut
            .loop_prescriptions
            .iter()
            .map(|value| value.parse())
            .collect::<Result<Vec<cuts::LoopPrescription>>>()?;
        let prepared = input.prepare_cut_projection(
            cut.index,
            &point,
            symbol(&request.epsilon, ns)?,
            &channel,
            prescriptions,
            &backend,
            &options,
            &context,
        )?;
        let conditions = prepared
            .nonzero_conditions()
            .iter()
            .map(Atom::to_canonical_string)
            .collect::<Vec<_>>();
        if let Some(samples) = cut.epsilon_samples {
            let samples = samples
                .iter()
                .map(|value| exact_rational(value, ns))
                .collect::<CliResult<Vec<_>>>()?;
            let values = prepared.evaluate_samples(&samples, &options, &context)?;
            return Ok(json!({
                "schema_version": 1, "operation": "cut_graph_samples", "cut_index": cut.index,
                "epsilon_samples": samples.iter().map(ToString::to_string).collect::<Vec<_>>(),
                "values": values.iter().map(|row| complex_json(&row[0])).collect::<Vec<_>>(),
                "working_bits": Precision::decimal(options.digits + options.guard_digits)?.bits,
                "verified_digits": Value::Null, "nonzero_conditions": conditions,
            }));
        }
        (
            prepared
                .solve(request.last_epsilon_power, &options, &context)?
                .remove(0),
            "cut_graph",
            Some(conditions),
        )
    } else {
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
        (result, "graph", None)
    };
    let mut output = json!({
        "schema_version":1,"operation":operation,"verified_digits":result.verified_digits,
        "working_bits":result.working_bits,"samples":result.samples,
        "validation_samples":result.validation_samples,"refinements":result.refinements,
        "coefficients":result.coefficients.iter().map(|(k,v)| (k.to_string(),complex_json(v))).collect::<BTreeMap<_,_>>(),
        "comparison_errors":result.comparison_errors.iter().map(|(k,v)| (k.to_string(),v.as_raw().to_string())).collect::<BTreeMap<_,_>>()
    });
    if let Some(conditions) = conditions {
        output["nonzero_conditions"] = json!(conditions);
    }
    Ok(output)
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
    /// Supplied original-basis range; defaults to the requested output last order.
    #[serde(default)]
    last_epsilon_power: Option<i32>,
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

#[derive(Clone, Copy, Deserialize)]
enum ContinuationSide {
    #[serde(rename = "+i0")]
    Plus,
    #[serde(rename = "-i0")]
    Minus,
}
impl From<ContinuationSide> for Prescription {
    fn from(side: ContinuationSide) -> Self {
        match side {
            ContinuationSide::Plus => Self::PlusI0,
            ContinuationSide::Minus => Self::MinusI0,
        }
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct PolynomialSideInput {
    polynomial: String,
    side: ContinuationSide,
}
#[derive(Deserialize)]
#[serde(rename_all = "snake_case")]
enum AdmissionInput {
    AllPlannerRoutesInDeclaredDomain,
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum ContinuationInput {
    PrescribedAffine {
        domain: String,
        prescriptions: Vec<PolynomialSideInput>,
        unprescribed_side: ContinuationSide,
        homotopy_admission: AdmissionInput,
    },
}
impl ContinuationInput {
    fn native(&self, namespace: &str) -> CliResult<PhysicalContinuation> {
        let Self::PrescribedAffine {
            domain,
            prescriptions,
            unprescribed_side,
            homotopy_admission: AdmissionInput::AllPlannerRoutesInDeclaredDomain,
        } = self;
        Ok(PhysicalContinuation {
            domain: domain.clone(),
            prescriptions: prescriptions
                .iter()
                .map(|p| {
                    Ok(PolynomialPrescription {
                        polynomial: parse(&p.polynomial, namespace)?,
                        prescription: p.side.into(),
                    })
                })
                .collect::<CliResult<_>>()?,
            unprescribed_side: (*unprescribed_side).into(),
        })
    }
}

// This policy is selected only by the explicit, mandatory declaration above.
// Local root signs and the contour planner cannot certify global monodromy.
fn admit_declared_route(_: &CachedBoundary, _: &CachedPoint, _: &PhysicalRoute) -> Result<bool> {
    Ok(true)
}
fn bind_continuation<S>(flow: RustFlow<S>, request: &TransportRequest) -> CliResult<RustFlow<S>> {
    match &request.continuation {
        Some(input) => Ok(flow.with_prescribed_continuation(input.native(&request.namespace)?)?),
        None => Ok(flow),
    }
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct TransportRequest {
    schema_version: u32,
    #[serde(default)]
    epsilon_shearing: bool,
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
    continuation: Option<ContinuationInput>,
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
    if request.epsilon_shearing && request.canonical.is_some() {
        return Err("epsilon_shearing accepts dense rational or registered-root derivatives; canonical inputs are already epsilon regular and are unsupported for explicit shearing".into());
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
        let flow = bind_continuation(flow, &request)?;
        let mut result = execute_transport(
            &request,
            directory,
            flow.identity(),
            range,
            None,
            |cache, point, policy| {
                if request.continuation.is_some() {
                    Ok(flow.evaluate_prescribed_to(
                        cache,
                        &point.restart_coordinates()?,
                        point.root_germ(),
                        range,
                        &options,
                        &context,
                        policy,
                        &admit_declared_route,
                    )?)
                } else {
                    Ok(flow.evaluate_to(
                        cache,
                        &point.restart_coordinates()?,
                        point.root_germ(),
                        range,
                        &options,
                        &context,
                        policy,
                    )?)
                }
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
        let flow = bind_continuation(flow, &request)?;
        if request.epsilon_shearing {
            let sheared = flow.regularize_epsilon(&context)?;
            let cached_range = sheared.required_range(range)?;
            return execute_transport(
                &request,
                directory,
                sheared.flow().identity(),
                range,
                Some(&sheared),
                |cache, point, policy| {
                    if request.continuation.is_some() {
                        Ok(sheared.flow().evaluate_prescribed_to(
                            cache,
                            &point.restart_coordinates()?,
                            cached_range,
                            &options,
                            &context,
                            policy,
                            &admit_declared_route,
                        )?)
                    } else {
                        Ok(sheared.flow().evaluate_to(
                            cache,
                            &point.restart_coordinates()?,
                            cached_range,
                            &options,
                            &context,
                            policy,
                        )?)
                    }
                },
            );
        }
        execute_transport(
            &request,
            directory,
            flow.identity(),
            range,
            None,
            |cache, point, policy| {
                if request.continuation.is_some() {
                    Ok(flow.evaluate_prescribed_to(
                        cache,
                        &point.restart_coordinates()?,
                        range,
                        &options,
                        &context,
                        policy,
                        &admit_declared_route,
                    )?)
                } else {
                    Ok(flow.evaluate_to(
                        cache,
                        &point.restart_coordinates()?,
                        range,
                        &options,
                        &context,
                        policy,
                    )?)
                }
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
        let flow = bind_continuation(flow, &request)?;
        if request.epsilon_shearing {
            let sheared = flow.regularize_epsilon(&context)?;
            let cached_range = sheared.required_range(range)?;
            return execute_transport(
                &request,
                directory,
                sheared.flow().identity(),
                range,
                Some(&sheared),
                |cache, point, policy| {
                    let germ = point.root_germ().ok_or("missing destination root_germ")?;
                    if request.continuation.is_some() {
                        Ok(sheared.flow().evaluate_prescribed_to(
                            cache,
                            &point.restart_coordinates()?,
                            germ,
                            cached_range,
                            &options,
                            &context,
                            policy,
                            &admit_declared_route,
                        )?)
                    } else {
                        Ok(sheared.flow().evaluate_to(
                            cache,
                            &point.restart_coordinates()?,
                            germ,
                            cached_range,
                            &options,
                            &context,
                            policy,
                        )?)
                    }
                },
            );
        }

        execute_transport(
            &request,
            directory,
            flow.identity(),
            range,
            None,
            |cache, point, policy| {
                if request.continuation.is_some() {
                    Ok(flow.evaluate_prescribed_to(
                        cache,
                        &point.restart_coordinates()?,
                        point.root_germ().ok_or("missing destination root_germ")?,
                        range,
                        &options,
                        &context,
                        policy,
                        &admit_declared_route,
                    )?)
                } else {
                    Ok(flow.evaluate_to(
                        cache,
                        &point.restart_coordinates()?,
                        point.root_germ().ok_or("missing destination root_germ")?,
                        range,
                        &options,
                        &context,
                        policy,
                    )?)
                }
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

// Type-erased boundary mapping only: all exact shifts and authentication remain
// in the shared library adapter, while evaluation stays in its typed flow.
trait CliEpsilonShearing {
    fn original_identity(&self) -> &BoundaryIdentity;
    fn shearing(&self) -> &symbolica_amflow::EpsilonShearing;
    fn required_range(&self, range: EpsilonRange) -> symbolica_amflow::Result<EpsilonRange>;
    fn to_sheared_boundary(
        &self,
        source: &CachedBoundary,
        last: i32,
    ) -> symbolica_amflow::Result<CachedBoundary>;
    fn to_original_boundary(
        &self,
        source: &CachedBoundary,
        range: EpsilonRange,
    ) -> symbolica_amflow::Result<CachedBoundary>;
}
impl<S> CliEpsilonShearing for EpsilonShearedFlow<S> {
    fn original_identity(&self) -> &BoundaryIdentity {
        EpsilonShearedFlow::original_identity(self)
    }
    fn shearing(&self) -> &symbolica_amflow::EpsilonShearing {
        EpsilonShearedFlow::shearing(self)
    }
    fn required_range(&self, range: EpsilonRange) -> symbolica_amflow::Result<EpsilonRange> {
        EpsilonShearedFlow::required_range(self, range)
    }
    fn to_sheared_boundary(
        &self,
        source: &CachedBoundary,
        last: i32,
    ) -> symbolica_amflow::Result<CachedBoundary> {
        EpsilonShearedFlow::to_sheared_boundary(self, source, last)
    }
    fn to_original_boundary(
        &self,
        source: &CachedBoundary,
        range: EpsilonRange,
    ) -> symbolica_amflow::Result<CachedBoundary> {
        EpsilonShearedFlow::to_original_boundary(self, source, range)
    }
}

// Both native connection kinds share parsing, bank I/O, progressive insertion,
// policy and result serialization. Numerical transport stays in RustFlow.
fn execute_transport(
    request: &TransportRequest,
    directory: &Path,
    identity: &BoundaryIdentity,
    range: EpsilonRange,
    sheared: Option<&dyn CliEpsilonShearing>,
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
    let source_identity = sheared.map_or(identity, CliEpsilonShearing::original_identity);
    for seed in &request.seeds {
        let p = Precision::decimal(seed.working_digits)?;
        let source = CachedBoundary {
            identity: source_identity.clone(),
            point: transport_point(
                &seed.coordinates,
                seed.root_germ.as_ref(),
                source_identity,
                ns,
            )?,
            kind: PointKind::Physical,
            range: EpsilonRange::new(range.leading, seed.last_epsilon_power.unwrap_or(range.last))?,
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
        };
        let source = if let Some(adapter) = sheared {
            let required = adapter.required_range(range)?;
            if source.range.last < required.last {
                return Err(format!(
                    "epsilon shearing requires original seed coefficients through epsilon^{}, but seed.last_epsilon_power is {}",
                    required.last, source.range.last
                ).into());
            }
            // Preserve the complete declared rectangular source range. Extra
            // supplied coefficients remain useful for future compatible queries.
            adapter.to_sheared_boundary(&source, source.range.last)?
        } else {
            source
        };
        seeds.push(source);
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
        let mut result = evaluate(&mut cache, &point, &policy)?;
        if let Some(adapter) = sheared {
            // The bank keeps its authenticated rescaled basis; only the result
            // is restored. Failure here leaves this destination unsaved.
            result.boundary = adapter.to_original_boundary(&result.boundary, range)?;
        }
        cache.save(&cache_directory)?;
        let mut output = json!({
            "coordinates":point_json(&result.boundary.point)?,
            "starting_point":point_json(&result.starting_point)?,
            "verified_digits":result.boundary.accuracy.verified_digits(),
            "working_bits":result.boundary.accuracy.working_bits(),
            "conditioning_digits":result.transport.as_ref().and_then(|t| t.diagnostics.conditioning_digits),
            "steps":result.transport.as_ref().map_or(0,|t| t.diagnostics.steps),
            "rejected_steps":result.transport.as_ref().map_or(0,|t| t.diagnostics.rejected_steps),
            "predicate_evaluations":result.transport.as_ref().map_or(0,|t| t.diagnostics.predicate_evaluations),
            "rational_trials":result.transport.as_ref().map_or(0,|t| t.diagnostics.pade_trials),
            "rational_steps":result.transport.as_ref().map_or(0,|t| t.diagnostics.pade_steps),
            "rational_fallbacks":result.transport.as_ref().map_or(0,|t| t.diagnostics.pade_fallbacks),
            "last_rational_fallback":result.transport.as_ref().and_then(|t| t.diagnostics.last_pade_fallback.clone()),
            "fundamental_boundary_charts":result.transport.as_ref().map_or(0,|t| t.diagnostics.fundamental_boundary_charts),
            "fundamental_boundary_fallback":result.transport.as_ref().and_then(|t| t.diagnostics.fundamental_boundary_fallback.clone()),
            "fundamental_boundary_retry":result.transport.as_ref().and_then(|t| t.diagnostics.fundamental_boundary_retry.clone()),
            "superseded_successes":result.transport.as_ref().map_or(0,|t| t.diagnostics.superseded_successes),
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
    let mut output = json!({"schema_version":1,"operation":"transport","loaded_boundaries":loaded,"retained_boundaries":cache.len(),"results":results});
    let output_identity = sheared.map_or(identity, CliEpsilonShearing::original_identity);
    if let Some(adapter) = sheared {
        let cached_range = adapter.required_range(range)?;
        output["identity"] = output_identity.key().into();
        output["epsilon_shearing"] = json!({
            "weights": adapter.shearing().weights(),
            "convention": "original_I_i = epsilon^weight_i * cached_J_i",
            "original_identity": adapter.original_identity().key(),
            "cached_identity": identity.key(),
            "original_output_range": {"leading":range.leading,"last":range.last},
            "requested_cached_range": {"leading":cached_range.leading,"last":cached_range.last},
        });
    }
    if let Some(continuation) = identity.physical_continuation() {
        let side = |s: Prescription| match s {
            Prescription::PlusI0 => "+i0",
            Prescription::MinusI0 => "-i0",
        };
        output["identity"] = output_identity.key().into();
        output["continuation"] = json!({
            "kind": "prescribed_affine",
            "domain": continuation.domain,
            "prescriptions": continuation.prescriptions.iter().map(|p| json!({
                "polynomial": p.polynomial.to_canonical_string(),
                "side": side(p.prescription),
            })).collect::<Vec<_>>(),
            "unprescribed_side": side(continuation.unprescribed_side),
            "homotopy_admission": "all_planner_routes_in_declared_domain",
        });
    }
    Ok(output)
}

fn run() -> CliResult<()> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if arguments.is_empty() || matches!(arguments[0].as_str(), "--help" | "-h" | "help") {
        println!(
            "RustFlow\n\nUsage: rustflow graph INPUT.json\n       rustflow transport INPUT.json\n       rustflow finite-density-prepare INPUT.json\n       rustflow finite-density INPUT.json\n       rustflow finite-density-sample INPUT.json EPSILON\n\nGraph input references native HEPKit model JSON and DOT files.\nTransport input supplies differential equations and a binary cache directory.\nFinite-density preparation validates graph/charge evidence and converts targets.\nFinite-density evaluation requires every cut sector to be admitted and closed; unsupported or unresolved sectors fail explicitly.\nFinite-density-sample uses an exact nonzero rational epsilon, for example 4/5.\nFinite-density results use the unscaled Euclidean measure; see docs/finite-density.md for current limits.\nPaths are relative to INPUT.json; results are JSON on stdout.\nSee docs/cli.md for exact-input and accuracy conventions."
        );
        return Ok(());
    }
    if arguments == ["--version"] {
        println!("RustFlow {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    let required_arguments = if arguments[0] == "finite-density-sample" {
        3
    } else {
        2
    };
    if arguments.len() != required_arguments {
        return Err("expected a command and one steering JSON file (plus exact epsilon for finite-density-sample); use --help".into());
    }
    let path = Path::new(&arguments[1]);
    let directory = path.parent().unwrap_or_else(|| Path::new("."));
    let text = fs::read_to_string(path)?;
    let output = match arguments[0].as_str() {
        "graph" => graph(serde_json::from_str(&text)?, directory)?,
        "transport" => transport(serde_json::from_str(&text)?, directory)?,
        "finite-density-prepare" => {
            let input: finite_density::DensityInput = serde_json::from_str(&text)?;
            input.prepare()?.preparation_report(65536)?
        }
        "finite-density" | "finite-density-sample" => {
            let input: finite_density::DensityInput = serde_json::from_str(&text)?;
            let epsilon = arguments.get(2).map(String::as_str);
            finite_density::interface::evaluate_request(
                &input,
                epsilon,
                &finite_density::interface::default_options(&input),
                finite_density::interface::default_closure_options(),
                &RunContext::default(),
                8,
            )?
        }
        _ => {
            return Err(
                "unknown command; expected graph, transport, finite-density-prepare, finite-density or finite-density-sample".into(),
            );
        }
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

#[cfg(test)]
mod step_size_settings_tests {
    use super::*;
    #[test]
    fn automatic_boundary_errors_are_default_and_validated() {
        assert_eq!(
            Settings::default()
                .options(Path::new("."))
                .unwrap()
                .boundary_error_strategy,
            BoundaryErrorStrategy::Automatic
        );
        let settings: Settings =
            serde_json::from_value(json!({"boundary_error_strategy":"fundamental_matrix"}))
                .unwrap();
        assert_eq!(
            settings
                .options(Path::new("."))
                .unwrap()
                .boundary_error_strategy,
            BoundaryErrorStrategy::FundamentalMatrix
        );
        let invalid: Settings =
            serde_json::from_value(json!({"boundary_error_strategy":"unchecked"})).unwrap();
        assert!(invalid.options(Path::new(".")).is_err());
    }
    #[test]
    fn integer_residual_arithmetic_is_explicit_and_validated() {
        assert_eq!(
            Settings::default()
                .options(Path::new("."))
                .unwrap()
                .residual_arithmetic,
            ResidualArithmetic::Ball
        );
        for (value, expected) in [
            ("ball", ResidualArithmetic::Ball),
            ("adaptive_integer", ResidualArithmetic::AdaptiveInteger),
        ] {
            let settings: Settings =
                serde_json::from_value(json!({"residual_arithmetic":value})).unwrap();
            assert_eq!(
                settings
                    .options(Path::new("."))
                    .unwrap()
                    .residual_arithmetic,
                expected
            );
        }
        let invalid: Settings =
            serde_json::from_value(json!({"residual_arithmetic":"integer"})).unwrap();
        assert!(invalid.options(Path::new(".")).is_err());
    }
    #[test]
    fn rational_trials_are_opt_in_and_degree_is_validated() {
        assert!(
            Settings::default()
                .options(Path::new("."))
                .unwrap()
                .pade
                .is_none()
        );
        let requested: Settings = serde_json::from_value(json!({"pade_degree":8})).unwrap();
        assert_eq!(
            requested
                .options(Path::new("."))
                .unwrap()
                .pade
                .unwrap()
                .degree,
            8
        );
        for degree in [0, 33] {
            let invalid: Settings = serde_json::from_value(json!({"pade_degree":degree})).unwrap();
            assert!(invalid.options(Path::new(".")).is_err());
        }
    }
    #[test]
    fn step_strategy_is_explicit_and_unknown_values_are_rejected() {
        let defaults = Settings::default().options(Path::new(".")).unwrap();
        assert_eq!(defaults.step_size_strategy, StepSizeStrategy::Halving);
        for (text, expected) in [
            ("halving", StepSizeStrategy::Halving),
            ("bracketed", StepSizeStrategy::Bracketed),
        ] {
            let settings: Settings =
                serde_json::from_value(json!({"step_size_strategy":text})).unwrap();
            assert_eq!(
                settings.options(Path::new(".")).unwrap().step_size_strategy,
                expected
            );
        }
        let invalid: Settings =
            serde_json::from_value(json!({"step_size_strategy":"unchecked"})).unwrap();
        assert!(invalid.options(Path::new(".")).is_err());
        assert_eq!(defaults.local_coordinate, LocalCoordinate::Identity);
        for (text, expected) in [
            ("identity", LocalCoordinate::Identity),
            ("balanced_mobius", LocalCoordinate::BalancedMobius),
        ] {
            let settings: Settings =
                serde_json::from_value(json!({"local_coordinate": text})).unwrap();
            assert_eq!(
                settings.options(Path::new(".")).unwrap().local_coordinate,
                expected
            );
        }
        let invalid: Settings =
            serde_json::from_value(json!({"local_coordinate": "unchecked"})).unwrap();
        assert!(invalid.options(Path::new(".")).is_err());
    }
}
