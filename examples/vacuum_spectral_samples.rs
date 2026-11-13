//! Epsilon samples of a vacuum family through the RustRed LiteRed lane, with
//! optional spectral gluing. Reproduces the four-loop P14 calibration and the
//! glued five-loop samples of the vacuum-master campaign from
//! `fixtures/vacuum/p14-family.json` (and a gluing map).
//!
//! ```text
//! vacuum_spectral_samples sample   --family F.json --dir OUT [--attempt A] [--indices 0,5] [--workers N] [options]
//! vacuum_spectral_samples spectral --family F.json --map MAP.json --dir OUT [--attempt A] [--indices ..] [options]
//! vacuum_spectral_samples fit      --dir OUT [--attempt A]
//! ```
//!
//! Options (defaults in brackets): `--digits` [20], `--guard` [116],
//! `--arithmetic-guard` [16], `--order` [176], `--deform-slot` [1],
//! `--max-steps` [20000], `--lane-workers` [2], `--child-fallback native|none`
//! [native], `--normalization-loops` [4; spectral: 5], `--arm-c` [5/2],
//! `--arm-t0` [120], `--arm-tcheck` [64], `--arm-digits` [working - 10].
//!
//! Attempt `A` uses RustFlow's `EpsilonFitProfile` for leading power -8: its
//! symmetric grid, guard digits raised to the profile's working digits and the
//! series order raised by 16 per attempt. Samples are written atomically as
//! `OUT/attempt-A/sample-III.json`; existing samples are skipped.
use serde::Deserialize;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
    time::Instant,
};
use symbolica::domains::float::Real;
use symbolica::prelude::*;
use symbolica_amflow::{
    epsilon::EpsilonFitProfile,
    recursive::RecursiveBoundary,
    spectral::{self, gluing},
    *,
};

const LEADING: i32 = -8;
const LAST: i32 = 0;
type AnyResult<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

#[derive(Deserialize)]
struct FamilyFile {
    name: String,
    epsilon: String,
    loops: Vec<String>,
    dimension: i64,
    physical_propagators: usize,
    denominators: Vec<Denominator>,
    targets: Vec<Vec<i16>>,
}

#[derive(Deserialize)]
struct Denominator {
    constant: String,
    scalar_products: Vec<String>,
}

fn atom(text: &str) -> AnyResult<Atom> {
    Ok(Atom::parse(
        text,
        "feynkit_graph",
        ParseSettings::default(),
    )?)
}

fn load_family(path: &Path) -> AnyResult<(IntegralFamily, Vec<Integral>)> {
    let file: FamilyFile = serde_json::from_slice(&fs::read(path)?)?;
    let propagators = file
        .denominators
        .iter()
        .map(|d| {
            Ok(Propagator {
                constant: atom(&d.constant)?,
                scalar_products: d
                    .scalar_products
                    .iter()
                    .map(|s| atom(s))
                    .collect::<AnyResult<Vec<_>>>()?,
            })
        })
        .collect::<AnyResult<Vec<_>>>()?;
    let family = IntegralFamily {
        name: file.name,
        loops: file.loops,
        external: vec![],
        external_gram: vec![],
        propagators,
        physical_propagators: file.physical_propagators,
        epsilon: symbol!(&file.epsilon),
        dimension: file.dimension,
    };
    let targets = file.targets.into_iter().map(Integral).collect();
    Ok((family, targets))
}

/// Digest of the family's denominators (as recorded by the campaign driver).
fn family_digest(f: &IntegralFamily) -> String {
    let text = serde_json::to_string(&json!({
        "loops": f.loops.len(),
        "physical": f.physical_propagators,
        "props": f.propagators.iter().map(|p| json!([p.constant.to_canonical_string(),
            p.scalar_products.iter().map(|a| a.to_canonical_string()).collect::<Vec<_>>()])).collect::<Vec<_>>(),
    }))
    .unwrap();
    blake3::hash(text.as_bytes()).to_hex()[..16].to_string()
}

struct Args(BTreeMap<String, String>);
impl Args {
    fn parse(items: impl Iterator<Item = String>) -> AnyResult<Self> {
        let mut map = BTreeMap::new();
        let mut items = items.peekable();
        while let Some(key) = items.next() {
            let key = key
                .strip_prefix("--")
                .ok_or_else(|| format!("unexpected argument {key}"))?
                .to_string();
            let value = match items.peek() {
                Some(v) if !v.starts_with("--") => items.next().unwrap(),
                _ => String::new(),
            };
            map.insert(key, value);
        }
        Ok(Self(map))
    }
    fn get(&self, key: &str) -> Option<&str> {
        self.0.get(key).map(String::as_str)
    }
    fn num<T: std::str::FromStr>(&self, key: &str, default: T) -> AnyResult<T>
    where
        T::Err: std::fmt::Display,
    {
        match self.get(key) {
            Some(v) => v
                .parse()
                .map_err(|e: T::Err| format!("--{key}: {e}").into()),
            None => Ok(default),
        }
    }
}

fn cjson(v: &ComplexFloat) -> Value {
    json!({"re": v.re.to_string(), "im": v.im.to_string()})
}

fn atomic(path: &Path, value: &Value) -> AnyResult<()> {
    let tmp = path.with_extension(format!("tmp-{}", std::process::id()));
    fs::write(&tmp, serde_json::to_vec_pretty(value)?)?;
    fs::rename(tmp, path)?;
    Ok(())
}

fn peak_rss_kib() -> u64 {
    fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1))
                .and_then(|v| v.parse().ok())
        })
        .unwrap_or(0)
}

fn options(a: &Args) -> AnyResult<FlowOptions> {
    let options = FlowOptions {
        digits: a.num("digits", 20)?,
        guard_digits: a.num("guard", 116)?,
        arithmetic_guard_digits: a.num("arithmetic-guard", 16)?,
        series_order: a.num("order", 176)?,
        workers: 1,
        sampled_reduction: true,
        mass_mode: MassMode::Propagators(vec![a.num("deform-slot", 1)?]),
        reuse_samples: false,
        max_steps: a.num("max-steps", 20000)?,
        max_boundary_attempts: 8,
        ..Default::default()
    };
    options.validate()?;
    Ok(options)
}

/// Profile grid and options of refinement attempt `attempt`.
fn attempt(base: &FlowOptions, attempt: usize) -> AnyResult<(EpsilonFitProfile, FlowOptions)> {
    let profile = EpsilonFitProfile::new(base.digits, base.guard_digits, LEADING, LAST, attempt)?;
    let mut refined = base.clone();
    refined.guard_digits = profile.working_digits - base.digits;
    refined.series_order += attempt * 16;
    refined.validate()?;
    Ok((profile, refined))
}

type Backend = FallbackBackend<LiteRedRustRedBackend, RustRedBackend>;

fn backend(a: &Args, family: &IntegralFamily) -> AnyResult<Backend> {
    let fallback = match a.get("child-fallback").unwrap_or("native") {
        "none" => None,
        // Native search with its defaults (depth 2, 4096 targets, one worker).
        "native" => Some(RustRedBackend::default()),
        other => return Err(format!("unknown --child-fallback {other}").into()),
    };
    Ok(FallbackBackend {
        primary: LiteRedRustRedBackend::new(LiteRedBackendOptions {
            workers: a.num("lane-workers", 2)?,
            ..Default::default()
        }),
        fallback,
        fallback_below_loops: family.loops.len(),
    })
}

fn basis_digest(prepared: &PreparedFlow) -> AnyResult<String> {
    let basis: Vec<&Vec<i16>> = prepared.reduced.basis.iter().map(|i| &i.0).collect();
    Ok(blake3::hash(&serde_json::to_vec(&basis)?).to_hex()[..16].to_string())
}

struct Job {
    index: usize,
    epsilon: Rational,
}

/// Run `jobs` on `workers` plain threads with large stacks (deep symbolic and
/// series recursion); each job writes its own file.
fn run_jobs(
    jobs: &[Job],
    workers: usize,
    run: impl Fn(&Job) -> AnyResult<()> + Sync,
) -> Vec<String> {
    let next = std::sync::atomic::AtomicUsize::new(0);
    let errors = std::sync::Mutex::new(Vec::new());
    std::thread::scope(|scope| {
        for w in 0..workers.clamp(1, jobs.len().max(1)) {
            std::thread::Builder::new()
                .name(format!("sample-worker-{w}"))
                .stack_size(256 << 20)
                .spawn_scoped(scope, || {
                    loop {
                        let k = next.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
                        let Some(job) = jobs.get(k) else { break };
                        if let Err(e) = run(job) {
                            let text = format!("sample {} (eps {}): {e}", job.index, job.epsilon);
                            eprintln!("{text}");
                            errors.lock().unwrap().push(text);
                        }
                    }
                })
                .expect("spawn sample worker");
        }
    });
    errors.into_inner().unwrap()
}

fn mode_samples(a: &Args, spectral_mode: bool) -> AnyResult<()> {
    let (family, mut targets) = load_family(Path::new(a.get("family").ok_or("--family")?))?;
    let map: Option<gluing::GluingMap> = if spectral_mode {
        let map: gluing::GluingMap =
            serde_json::from_slice(&fs::read(a.get("map").ok_or("--map")?)?)?;
        targets = map.targets();
        Some(map)
    } else {
        None
    };
    let base = options(a)?;
    let attempt_index = a.num("attempt", 0usize)?;
    let (profile, refined) = attempt(&base, attempt_index)?;
    let dir = PathBuf::from(a.get("dir").ok_or("--dir")?).join(format!("attempt-{attempt_index}"));
    fs::create_dir_all(&dir)?;
    let working = refined.digits + refined.guard_digits;
    let deformed_slot: usize = a.num("deform-slot", 1)?;
    let arm_spec = gluing::ArmSpec {
        c: a.get("arm-c").unwrap_or("5/2").into(),
        t0: a.num("arm-t0", 120)?,
        t_check: a.num("arm-tcheck", 64)?,
        digits: a.num("arm-digits", working.saturating_sub(10))?,
    };
    let normalization_loops: i64 =
        a.num("normalization-loops", if spectral_mode { 5 } else { 4 })?;
    atomic(
        &dir.join("profile.json"),
        &json!({"attempt": attempt_index, "samples": profile.samples.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "radius": profile.radius.to_string(), "working_digits": working, "options": format!("{refined:?}"),
            "family_digest": family_digest(&family), "targets": targets.iter().map(|t| &t.0).collect::<Vec<_>>(),
            "arm": spectral_mode.then(|| json!({"c": arm_spec.c, "t0": arm_spec.t0, "t_check": arm_spec.t_check,
                "digits": arm_spec.digits}))}),
    )?;
    let wanted: Option<Vec<usize>> = a
        .get("indices")
        .map(|list| list.split(',').map(str::parse).collect())
        .transpose()?;
    let jobs: Vec<Job> = profile
        .samples
        .iter()
        .cloned()
        .enumerate()
        .filter(|(i, _)| wanted.as_ref().is_none_or(|w| w.contains(i)))
        .filter(|(i, _)| !dir.join(format!("sample-{i:03}.json")).exists())
        .map(|(index, epsilon)| Job { index, epsilon })
        .collect();
    let backend = backend(a, &family)?;
    let context = RunContext::default();
    let errors = run_jobs(&jobs, a.num("workers", 1)?, |job| {
        let t0 = Instant::now();
        let prepared = PreparedFlow::new_at_epsilon(
            &family,
            &targets,
            &KinematicPoint::default(),
            &backend,
            &refined,
            &context,
            &job.epsilon,
        )?;
        let prepare_seconds = t0.elapsed().as_secs_f64();
        let boundary = RecursiveBoundary::new(&backend, &refined, &context);
        let p = Precision::decimal(working)?;
        let mut out = json!({"epsilon": job.epsilon.to_string(), "index": job.index, "attempt": attempt_index,
            "working_digits": working, "working_bits": p.bits, "series_order": refined.series_order,
            "basis_size": prepared.reduced.basis.len(), "basis_digest": basis_digest(&prepared)?,
            "prepare_seconds": prepare_seconds});
        let t1 = Instant::now();
        if let Some(map) = &map {
            let ctx = spectral::Ctx { bits: p.bits };
            let arm = gluing::SpectralArm::new(&ctx, arm_spec.clone())?;
            let sample = prepared.evaluate_spectral(
                &job.epsilon,
                &refined,
                &boundary,
                &context,
                &arm.eta_path(&ctx),
                &arm.eta_nodes(),
                refined.series_order,
            )?;
            out["evaluate_seconds"] = json!(t1.elapsed().as_secs_f64());
            let t2 = Instant::now();
            let glued = gluing::glue(
                &ctx,
                &job.epsilon.to_string(),
                &arm,
                refined.series_order,
                &targets,
                deformed_slot,
                &sample,
                map,
                normalization_loops,
            )?;
            out["assembly_seconds"] = json!(t2.elapsed().as_secs_f64());
            out["endpoint"] = json!(
                targets
                    .iter()
                    .zip(&sample.endpoint)
                    .map(|(t, v)| json!({"target": t.0, "value": cjson(v)}))
                    .collect::<Vec<_>>()
            );
            out["masters"] = json!(glued.masters.iter().map(|(name, m)| (name.clone(), json!({"raw": spectral::fmt(&m.raw),
                "normalized": {"re": spectral::fmt(&m.normalized), "im": "0"}, "tail_split_digits": m.tail_split_digits})))
                .collect::<serde_json::Map<_, _>>());
            out["X"] = json!(glued.contour_integrals.iter().map(|x| json!({"key": x.key, "X": spectral::fmt(&x.x),
                "X_tail_split_digits": x.tail_split_digits,
                "integer_z_inverse_coefficient": x.integer_z_inverse_coefficient.as_ref().map(|c| json!({
                    "re": spectral::fmt_digits(c.real(), 20), "im": spectral::fmt_digits(c.imag(), 20)}))})).collect::<Vec<_>>());
            out["target_diagnostics"] = json!(glued.diagnostics.iter().map(|d| json!({"target": d.target,
                "classes": d.classes, "asym_vs_node_at_T0_digits": d.asym_vs_node_at_t0_digits,
                "eta_independent_node_vs_endpoint_digits": d.eta_independent_node_vs_endpoint_digits})).collect::<Vec<_>>());
            out["bubble_self_check_digits"] = json!(glued.bubble_self_check_digits);
            out["log_terms_max"] = json!(sample.log_terms_max);
            out["arm_nodes"] = json!(arm.arm.z.len());
            out["arm_panels"] = json!(format!("{:?}", arm.arm.panels));
            out["arm_segments"] = json!(sample.arm_segments);
            out["arm_steps"] = json!(sample.arm_steps);
            out["arm_rejected_steps"] = json!(sample.arm_rejected_steps);
        } else {
            let values = prepared.evaluate(&job.epsilon, &refined, &boundary, &context)?;
            out["evaluate_seconds"] = json!(t1.elapsed().as_secs_f64());
            let euler = ComplexFloat::new(p.real(0).euler(), p.real(0));
            let norm = p.exp(&p.mul(
                &p.rational(&(Rational::from(normalization_loops) * &job.epsilon)),
                &euler,
            ));
            out["endpoint"] = json!(targets.iter().zip(&values)
                .map(|(t, v)| json!({"target": t.0, "value": cjson(v), "normalized": cjson(&p.mul(v, &norm))}))
                .collect::<Vec<_>>());
        }
        out["process_peak_rss_kib"] = json!(peak_rss_kib());
        atomic(&dir.join(format!("sample-{:03}.json", job.index)), &out)?;
        Ok(())
    });
    if !errors.is_empty() {
        return Err(format!("{} samples failed; first: {}", errors.len(), errors[0]).into());
    }
    Ok(())
}

/// Laurent fit of every endpoint target's normalized value over the attempt grid.
fn mode_fit(a: &Args) -> AnyResult<()> {
    let attempt_index = a.num("attempt", 0usize)?;
    let dir = PathBuf::from(a.get("dir").ok_or("--dir")?).join(format!("attempt-{attempt_index}"));
    let profile: Value = serde_json::from_slice(&fs::read(dir.join("profile.json"))?)?;
    let working = profile["working_digits"].as_u64().ok_or("working_digits")? as u32;
    let p = Precision::decimal(working)?;
    let mut eps = vec![];
    let mut values: BTreeMap<String, Vec<ComplexFloat>> = BTreeMap::new();
    for (i, e) in profile["samples"]
        .as_array()
        .ok_or("samples")?
        .iter()
        .enumerate()
    {
        let sample: Value =
            serde_json::from_slice(&fs::read(dir.join(format!("sample-{i:03}.json")))?)?;
        let text = e.as_str().ok_or("epsilon")?;
        if sample["epsilon"] != text {
            return Err(format!("sample {i}: epsilon mismatch").into());
        }
        let (n, d) = text.split_once('/').unwrap_or((text, "1"));
        eps.push(Rational::from((n.parse::<i64>()?, d.parse::<i64>()?)));
        for entry in sample["endpoint"].as_array().ok_or("endpoint")? {
            let value = &entry["normalized"];
            values
                .entry(entry["target"].to_string())
                .or_default()
                .push(p.parse(
                    value["re"].as_str().ok_or("re")?,
                    value["im"].as_str().ok_or("im")?,
                )?);
        }
    }
    let mut out = serde_json::Map::new();
    for (target, v) in values {
        let fit = fit_epsilon(&eps, &v, LEADING, LAST, p)?;
        out.insert(
            target,
            json!(
                fit.coefficients
                    .iter()
                    .map(|(k, c)| (k.to_string(), cjson(c)))
                    .collect::<serde_json::Map<_, _>>()
            ),
        );
    }
    let result = json!({"attempt": attempt_index, "nodes": eps.len(), "working_digits": working, "coefficients": out});
    atomic(&dir.join("fit.json"), &result)?;
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}

fn main() -> AnyResult<()> {
    let mut argv = std::env::args().skip(1);
    let mode = argv.next().unwrap_or_default();
    let a = Args::parse(argv)?;
    match mode.as_str() {
        "sample" => mode_samples(&a, false),
        "spectral" => mode_samples(&a, true),
        "fit" => mode_fit(&a),
        _ => Err("usage: vacuum_spectral_samples sample|spectral|fit --dir OUT ...".into()),
    }
}
