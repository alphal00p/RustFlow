//! Long native ggHg boundary runs. References are loaded only after native fitting.
//! Run with --help; no Python or Mathematica runtime is required.
#[path = "../tests/support/gg_hg_anchors.rs"]
mod anchors;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    sync::{
        Arc, Mutex,
        atomic::{AtomicUsize, Ordering},
        mpsc,
    },
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};
use symbolica_amflow::gg_hg::{HiggsJetIntegralSystem, PluginFamilyKind};
use symbolica_amflow::symbolica::prelude::*;
use symbolica_amflow::transport_cache::{BoundaryCache, CachedBoundary, RootGerm, RootSheet};
use symbolica_amflow::{
    Error, FlowOptions, Precision, Progress, Result, RunContext, RustRedBackend,
};

const NAMESPACE: &str = "rustflow_gg_hg_boundaries";
const HELP: &str = "Native verified ggHg canonical boundary generation
  cargo run --release --example gg_hg_boundaries -- [OPTIONS]

  --family pl|np               Topology (default pl)
  --configuration LABEL        Use a published physical SOURCE configuration
  --list-configurations        List this topology's eight physical source labels
  --describe                   Print exact input/options without numerical work
  --digits N                   Requested mixed coefficient accuracy (default 20)
  --guard-digits N              Additional working digits (default 40)
  --order N                    Initial series order (default 80)
  --workers N                  Independent epsilon workers (default 1)
  --last-epsilon N              Last fitted canonical power (default 4)
  --max-precision-attempts N    Refinement limit (default 3)
  --max-steps N                Continuation acceptance budget (default 1000)
  --case-batch N               Native reducer batch size (default 16)
  --cache-directory PATH        Bank, samples, exact cache, reports
                               (default target/gg-hg-boundaries)
  --checkpoint-directory PATH  Override native IBP checkpoint directory
  --force                      Recompute numerical samples and boundary
  --cancel-file PATH           Cooperatively cancel when PATH exists
  --help                       Show this message

Defaults are exact Euclidean PL (-1/10,-1/25,-1/50) or NP (-1/10,-1/5,-1),
all principal roots, D=4-2 epsilon, and +i0. Euclidean anchor coefficients are
compared ONLY AFTER native generation succeeds. No reference seed fallback.
One process per cache directory; --workers bounds sample concurrency within it.
";

struct Arguments {
    family: PluginFamilyKind,
    configuration: Option<String>,
    directory: PathBuf,
    checkpoints: Option<PathBuf>,
    cancel_file: Option<PathBuf>,
    options: FlowOptions,
    last: i32,
    batch: usize,
    force: bool,
    describe: bool,
    list: bool,
    help: bool,
}
impl Arguments {
    fn parse(arguments: impl IntoIterator<Item = String>) -> Result<Self> {
        let mut out = Self {
            family: PluginFamilyKind::Planar,
            configuration: None,
            directory: "target/gg-hg-boundaries".into(),
            checkpoints: None,
            cancel_file: None,
            options: FlowOptions::default(),
            last: 4,
            batch: 16,
            force: false,
            describe: false,
            list: false,
            help: false,
        };
        let mut arguments = arguments.into_iter();
        while let Some(flag) = arguments.next() {
            if matches!(
                flag.as_str(),
                "--force" | "--describe" | "--list-configurations" | "--help"
            ) {
                match flag.as_str() {
                    "--force" => out.force = true,
                    "--describe" => out.describe = true,
                    "--list-configurations" => out.list = true,
                    _ => out.help = true,
                }
                continue;
            }
            let value = arguments
                .next()
                .ok_or_else(|| Error::InvalidInput(format!("missing value for {flag}")))?;
            let integer = || {
                value.parse::<usize>().map_err(|_| {
                    Error::InvalidInput(format!("{flag} expects a nonnegative integer"))
                })
            };
            let u32_value = || {
                u32::try_from(integer()?)
                    .map_err(|_| Error::InvalidInput(format!("{flag} exceeds u32")))
            };
            match flag.as_str() {
                "--family" => {
                    out.family = match value.as_str() {
                        "pl" => PluginFamilyKind::Planar,
                        "np" => PluginFamilyKind::Nonplanar,
                        _ => return Err(Error::InvalidInput("--family expects pl or np".into())),
                    }
                }
                "--configuration" => out.configuration = Some(value),
                "--cache-directory" => out.directory = value.into(),
                "--checkpoint-directory" => out.checkpoints = Some(value.into()),
                "--cancel-file" => out.cancel_file = Some(value.into()),
                "--digits" => out.options.digits = u32_value()?,
                "--guard-digits" => out.options.guard_digits = u32_value()?,
                "--order" => out.options.series_order = integer()?,
                "--workers" => out.options.workers = integer()?,
                "--max-steps" => out.options.max_steps = integer()?,
                "--max-precision-attempts" => out.options.max_precision_attempts = integer()?,
                "--case-batch" => out.batch = integer()?,
                "--last-epsilon" => {
                    out.last = i32::try_from(integer()?)
                        .map_err(|_| Error::InvalidInput("epsilon order exceeds i32".into()))?
                }
                _ => return Err(Error::InvalidInput(format!("unknown option {flag}"))),
            }
        }
        out.options.validate()?;
        if out.batch == 0 || out.last > 32 {
            return Err(Error::InvalidInput(
                "positive case batch and epsilon order <=32 required".into(),
            ));
        }
        Ok(out)
    }
}

fn seconds(duration: Duration) -> String {
    format!("{}.{:09}", duration.as_secs(), duration.subsec_nanos())
}
fn write_json(path: &Path, data: &Value) -> Result<()> {
    let temporary = path.with_extension(format!("{}.tmp", std::process::id()));
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)?;
    serde_json::to_writer_pretty(&mut file, data).map_err(|e| Error::Cache(e.to_string()))?;
    file.write_all(b"\n")?;
    file.sync_all()?;
    fs::rename(temporary, path)?;
    Ok(())
}
fn exact(text: &str) -> Result<Atom> {
    Atom::parse(text, NAMESPACE, Default::default()).map_err(|e| Error::InvalidInput(e.to_string()))
}
fn point(
    input: &HiggsJetIntegralSystem,
    args: &Arguments,
) -> Result<(String, BTreeMap<Symbol, Atom>, RootGerm)> {
    if let Some(label) = &args.configuration {
        let case = input
            .configurations(NAMESPACE)?
            .into_iter()
            .find(|c| &c.label == label)
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "unknown configuration for {}: {label}",
                    args.family.id()
                ))
            })?;
        return Ok((case.label, case.start, case.germ));
    }
    // These are exact mathematical coordinates, independent of reference files.
    let values = match args.family {
        PluginFamilyKind::Planar => ["-1/10", "-1/25", "-1/50"],
        PluginFamilyKind::Nonplanar => ["-1/10", "-1/5", "-1"],
    };
    let coordinates = input
        .basis_map()
        .coordinates
        .into_iter()
        .zip(values)
        .map(|(s, v)| Ok((s, exact(v)?)))
        .collect::<Result<_>>()?;
    let germ = RootGerm {
        sheets: input
            .basis_map()
            .roots
            .iter()
            .map(|r| (r.symbol, RootSheet::Principal))
            .collect(),
    };
    Ok((format!("{}-euclidean", args.family.id()), coordinates, germ))
}
fn boundary_json(boundary: &CachedBoundary) -> Value {
    json!({"leading_epsilon_power":boundary.range.leading,"last_epsilon_power":boundary.range.last,
        "canonical_basis_indices":(1..=boundary.identity.dimension()).collect::<Vec<_>>(),
        "coefficients":boundary.coefficients.iter().map(|row|row.iter().map(|v|json!({"real":v.re.to_string(),"imaginary":v.im.to_string()})).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "absolute_errors":boundary.accuracy.comparison_errors().iter().map(|row|row.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "verified_digits":boundary.accuracy.verified_digits(),"working_bits":boundary.accuracy.working_bits(),
        "accuracy_metric":"10^(-digits)*max(1,abs(coefficient)); empirical refinement evidence, not interval bounds",
        "provenance":boundary.accuracy.provenance(),"identity":boundary.identity.key()})
}
fn compare_anchor(input: &HiggsJetIntegralSystem, boundary: &CachedBoundary) -> Result<Value> {
    // Deliberately called only AFTER native generation returns successfully.
    let p = Precision {
        bits: boundary
            .accuracy
            .working_bits()
            .max(Precision::decimal(80)?.bits),
    };
    let anchor = anchors::EuclideanAnchor::load(input, p)?;
    let expected = anchor
        .point
        .0
        .iter()
        .map(|(a, v)| match a.as_view() {
            AtomView::Var(s) => Ok((s.get_symbol(), v.clone())),
            _ => Err(Error::InvalidInput("non-symbolic anchor coordinate".into())),
        })
        .collect::<Result<BTreeMap<_, _>>>()?;
    if boundary.point.restart_coordinates()? != expected
        || boundary.point.root_germ() != Some(&anchor.root_germ)
    {
        return Err(Error::InvalidInput(
            "native result is not at the comparison anchor and germ".into(),
        ));
    }
    let mut maximum = p.real(0);
    let mut maximum_ratio = p.real(0);
    let mut failures = vec![];
    let mut checked = 0;
    let last = boundary.range.last.min(4);
    for power in 0..=last {
        for index in 0..input.basis_map().kind.dimension() {
            let row = (power - boundary.range.leading) as usize;
            let difference = p.norm(&p.sub(
                &p.round(&boundary.coefficients[row][index]),
                &anchor.coefficients[power as usize][index],
            ));
            let allowance = boundary.accuracy.comparison_errors()[row][index].clone()
                + anchor.absolute_errors[power as usize][index].clone();
            let ratio = difference.clone() / allowance.clone();
            if difference > maximum {
                maximum = difference.clone();
            }
            if ratio > maximum_ratio {
                maximum_ratio = ratio;
            }
            if difference > allowance {
                failures.push(json!({"epsilon_power":power,"canonical_index":index+1,"absolute_difference":difference.to_string(),"combined_allowance":allowance.to_string()}));
            }
            checked += 1;
        }
    }
    Ok(
        json!({"passed":failures.is_empty(),"checked_coefficients":checked,"last_reference_power":last,
        "reference_accuracy_cap_digits":anchor.verified_digits,"reference_provenance":anchor.provenance,
        "metric":"absolute complex difference <= native allowance + recorded reference allowance",
        "maximum_absolute_difference":maximum.to_string(),"maximum_error_ratio":maximum_ratio.to_string(),"failures":failures,
        "reference_loaded_after_native_success":true}),
    )
}
fn run(args: Arguments) -> Result<()> {
    if args.help {
        println!("{HELP}");
        return Ok(());
    }
    let total = Instant::now();
    let input = HiggsJetIntegralSystem::load(args.family, NAMESPACE)?;
    if args.list {
        for case in input.configurations(NAMESPACE)? {
            println!("{}", case.label);
        }
        return Ok(());
    }
    let (label, coordinates, germ) = point(&input, &args)?;
    let mut options = args.options.clone();
    options.cache_directory = Some(args.directory.join("exact"));
    options.sample_cache_directory = Some(args.directory.join("samples"));
    options.reuse_samples = !args.force;
    let checkpoints = args
        .checkpoints
        .clone()
        .unwrap_or_else(|| args.directory.join("ibp"));
    let mut report = json!({"schema":"rustflow-native-gg-hg-boundary-run-v1","status":"prepared","label":label,
        "family":args.family.id(),"namespace":NAMESPACE,
        "coordinates":coordinates.iter().map(|(&s,v)|(Atom::var(s).to_canonical_string(),v.to_canonical_string())).collect::<BTreeMap<_,_>>(),
        "root_germs":germ.sheets.iter().map(|(&s,sheet)|(Atom::var(s).to_canonical_string(),format!("{sheet:?}"))).collect::<BTreeMap<_,_>>(),
        "options":{"digits":options.digits,"guard_digits":options.guard_digits,"series_order":options.series_order,
            "workers":options.workers,"last_epsilon_power":args.last,"force_numerical_recomputation":args.force,
            "max_precision_attempts":options.max_precision_attempts,"max_steps":options.max_steps,"reducer_case_batch":args.batch},
        "cache_directory":args.directory,"sample_directory":options.sample_cache_directory,"ibp_checkpoints":checkpoints,
        "build":{"profile":if cfg!(debug_assertions){"dev"}else{"release"},"rustflow_source_digest":env!("PORT_SOURCE_DIGEST"),"dependency_source_digest":env!("DEPENDENCY_SOURCE_DIGEST")},
        "mathematical_source":input.basis_map().paper.url,"physical_map_source_sha256":input.basis_map().evidence.source_sha256,
        "canonical_identity":input.transport().identity().key(),"reference_loaded":false,"verified_bank_saved":false});
    if args.describe {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).map_err(|e| Error::Cache(e.to_string()))?
        );
        return Ok(());
    }
    let identifier = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|e| Error::InvalidInput(e.to_string()))?
        .as_nanos();
    let run_directory = args
        .directory
        .join("runs")
        .join(format!("{label}-{identifier}-{}", std::process::id()));
    fs::create_dir_all(&run_directory)?;
    let report_file = run_directory.join("result.json");
    write_json(&report_file, &report)?;
    let events = Arc::new(Mutex::new(File::create(
        run_directory.join("progress.jsonl"),
    )?));
    let progress_failure = Arc::new(Mutex::new(None::<String>));
    let sample_events = Arc::new(AtomicUsize::new(0));
    let mut context = RunContext::default();
    let token = context.cancellation.clone();
    let counters = sample_events.clone();
    let failed = progress_failure.clone();
    context.progress = Some(Arc::new(move |event| {
        if matches!(&event, Progress::Sample { .. }) {
            counters.fetch_add(1, Ordering::Relaxed);
        }
        if !matches!(&event, Progress::Step { .. }) {
            eprintln!("[{}s] {event:?}", seconds(total.elapsed()));
        }
        let payload =
            json!({"elapsed_seconds":seconds(total.elapsed()),"event":format!("{event:?}")});
        let written = events
            .lock()
            .map_err(|_| "progress lock poisoned".to_owned())
            .and_then(|mut file| writeln!(file, "{payload}").map_err(|e| e.to_string()));
        if let Err(message) = written {
            if let Ok(mut failure) = failed.lock() {
                *failure = Some(message);
            }
            token.cancel();
        }
    }));
    let (stop, stopped) = mpsc::channel::<()>();
    let monitor = if let Some(path) = args.cancel_file.clone() {
        if path.exists() {
            context.cancellation.cancel();
        }
        let token = context.cancellation.clone();
        Some(
            std::thread::Builder::new()
                .name("gg-hg-cancel".into())
                .spawn(move || {
                    loop {
                        if path.exists() {
                            token.cancel();
                            break;
                        }
                        if !matches!(
                            stopped.recv_timeout(Duration::from_millis(250)),
                            Err(mpsc::RecvTimeoutError::Timeout)
                        ) {
                            break;
                        }
                    }
                })?,
        )
    } else {
        None
    };
    let generation = Instant::now();
    let result = (|| -> Result<()> {
        let mut bank = BoundaryCache::load(&args.directory.join("boundaries"))?;
        report["initial_bank_entries"] = json!(bank.len());
        let backend = RustRedBackend {
            checkpoints: Some(checkpoints),
            max_sector_batch: args.batch,
            ..Default::default()
        };
        let boundary = input.generate_boundary(
            &mut bank,
            &coordinates,
            &germ,
            args.last,
            &options,
            &backend,
            &context,
            !args.force,
        )?;
        report["native_generation_seconds"] = json!(seconds(generation.elapsed()));
        boundary.validate()?;
        report["boundary"] = boundary_json(&boundary);
        if args.configuration.is_none() {
            let comparison = compare_anchor(&input, &boundary)?;
            let passed = comparison["passed"] == true;
            report["reference_loaded"] = json!(true);
            report["reference_comparison"] = comparison;
            if !passed {
                return Err(Error::Accuracy("native Laurent coefficients disagree with the independent Euclidean anchor allowances".into()));
            }
        } else {
            report["reference_comparison"] = json!({"status":"not_run","reason":"physical source configuration; no Euclidean reference applies"});
        }
        context.cancellation.check()?;
        bank.save(&args.directory.join("boundaries"))?;
        report["verified_bank_saved"] = json!(true);
        report["final_bank_entries"] = json!(bank.len());
        Ok(())
    })();
    drop(stop);
    if let Some(monitor) = monitor {
        monitor
            .join()
            .map_err(|_| Error::Io(std::io::Error::other("cancel monitor panicked")))?;
    }
    let result = if let Some(message) = progress_failure
        .lock()
        .map_err(|_| Error::Io(std::io::Error::other("progress error lock poisoned")))?
        .take()
    {
        Err(Error::Io(std::io::Error::other(message)))
    } else {
        result
    };
    report["elapsed_seconds"] = json!(seconds(total.elapsed()));
    report["sample_progress_events"] = json!(sample_events.load(Ordering::Relaxed));
    match &result {
        Ok(()) => report["status"] = json!("success"),
        Err(error) => {
            report["status"] = json!(if matches!(error, Error::Cancelled) {
                "cancelled"
            } else {
                "failed"
            });
            report["error"] = json!({"typed":format!("{error:?}"),"message":error.to_string()});
        }
    }
    write_json(&report_file, &report)?;
    println!("{}", report_file.display());
    result
}
fn main() {
    if let Err(error) = Arguments::parse(std::env::args().skip(1)).and_then(run) {
        eprintln!("{error}");
        std::process::exit(if matches!(error, Error::Cancelled) {
            130
        } else {
            1
        });
    }
}
