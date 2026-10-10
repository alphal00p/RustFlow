//! Validation-only ordinary recursive hard coefficient. No reference input.
use serde_json::json;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};
use symbolica::prelude::*;
use symbolica_amflow::*;

fn main() {
    let directory = std::path::PathBuf::from(std::env::args().nth(1).expect("output directory"));
    std::fs::create_dir_all(&directory).unwrap();
    let family = IntegralFamily {
        name: "independent_equal_mass_hard_boundary_gate".into(),
        loops: vec!["K".into(), "L".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![
            Propagator::quadratic(&[1, 0], &[], Atom::one(), &[]).unwrap(),
            Propagator::quadratic(&[0, 1], &[], Atom::one(), &[]).unwrap(),
            Propagator::quadratic(&[1, -1], &[], Atom::one(), &[]).unwrap(),
        ],
        physical_propagators: 3,
        epsilon: symbol!("hard_boundary_gate::epsilon"),
        dimension: 4,
    };
    let epsilon = Rational::from((4, 5));
    let target = Integral(vec![1, 1, 1]);
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    for (digits, order) in [(18, 60), (28, 60), (28, 80)] {
        let path = directory.join(format!("prediction-{digits}-{order}.json"));
        assert!(!path.exists(), "preserve existing native predictions");
        let preparations = Arc::new(AtomicUsize::new(0));
        let observed = preparations.clone();
        let context = RunContext {
            progress: Some(Arc::new(move |event| {
                if matches!(event, Progress::Prepared { .. }) {
                    observed.fetch_add(1, Ordering::Relaxed);
                }
                if !matches!(event, Progress::Step { .. }) {
                    eprintln!("{event:?}");
                }
            })),
            ..Default::default()
        };
        let options = FlowOptions {
            digits,
            guard_digits: 40,
            series_order: order,
            ..Default::default()
        };
        let precision = Precision::decimal(digits + options.guard_digits).unwrap();
        let provider = recursive::RecursiveBoundary::new(&backend, &options, &context)
            .with_terminal_policy(recursive::RecursiveTerminalPolicy::TadpolesOnly);
        let began = std::time::Instant::now();
        let value = match provider.evaluate(&family, &target, &epsilon, precision) {
            Ok(value) => value,
            Err(error) => {
                std::fs::write(
                    directory.join("failure.json"),
                    serde_json::to_vec_pretty(&json!({
                        "digits":digits,"series_order":order,"error":error.to_string(),
                        "elapsed_seconds":began.elapsed().as_secs_f64(),"numerical_acceptance":false
                    }))
                    .unwrap(),
                )
                .unwrap();
                panic!("recursive hard boundary failed: {error}");
            }
        };
        assert!(
            preparations.load(Ordering::Relaxed) > 0,
            "the nonterminal two-loop coefficient must run recursive AMF"
        );
        std::fs::write(
            path,
            serde_json::to_vec_pretty(&json!({
                "scope":"Ordinary recursive equal-mass hard coefficient; no reference data loaded",
                "family":format!("{family:?}"),"powers":target.0,"dimension":"12/5","epsilon":"4/5",
                "normalization":"d^D K/(i*pi^(D/2)) d^D L/(i*pi^(D/2))",
                "bubble_subloops":false,"terminal_policy":"TadpolesOnly",
                "digits":digits,"guard_digits":40,"series_order":order,
                "working_bits":precision.bits,"value":value.to_string(),
                "native_preparations":preparations.load(Ordering::Relaxed),
                "elapsed_seconds":began.elapsed().as_secs_f64(),"reference_files_read":0
            }))
            .unwrap(),
        )
        .unwrap();
    }
}
