//! Validation only: compare identical source discovery across native revisions.
use rustred::algebra::CoefficientContext;
use rustred::persistence::equivalent_generated_programs;
use rustred::solver::guarded::{GuardedProgram, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole};
use rustred::solver::{Integral, SearchOptions, Term};
use std::{path::Path, sync::Arc};

fn source() -> Arc<GuardedSourceSystem<1>> {
    let context = CoefficientContext::try_new(["priority_n", "priority_x"]).unwrap();
    let x = context.parameter("priority_x").unwrap();
    Arc::new(GuardedSourceSystem::new(
        "priority-domain-regression",
        [IndexRole::Ordinary],
        [0],
        vec![GuardedSource::new(
            "positive-ray-recurrence",
            vec![
                Term { integral: Integral::symbolic([0]).unwrap(), coefficient: context.one().numerator },
                Term { integral: Integral::symbolic([-1]).unwrap(), coefficient: (-x.clone()).numerator },
            ],
            IndexDomain::new([IndexBounds::new(Some(2), None).unwrap()]).unwrap(),
        ).with_nonzero_conditions(vec![x.numerator])],
    ).unwrap())
}

fn main() {
    let args = std::env::args().collect::<Vec<_>>();
    let output = Path::new(&args[1]);
    std::fs::create_dir_all(output).unwrap();
    for budget in [3, 32] {
        let sources = source();
        let found = sources.solve_domains(
            vec![IndexDomain::unrestricted()],
            SearchOptions { max_depth: Some(2), sample_seed: 0, ..Default::default() },
            budget,
        ).unwrap();
        let gaps = found.unresolved.iter().map(|gap|
            format!("{:?}\t{:?}\t{}\n", gap.domain, gap.reason, gap.detail)
        ).collect::<String>();
        let program = GuardedProgram::new(sources.clone(), found.rules, []).unwrap();
        let bytes = program.encode_native(Default::default()).unwrap();
        // Both old and new executions independently replay all sealed rules.
        let decoded = GuardedProgram::decode_generated(&bytes, sources, Default::default()).unwrap();
        let application = decoded.apply(&[2]).unwrap();
        std::fs::write(output.join(format!("budget-{budget}.bin")), &bytes).unwrap();
        std::fs::write(output.join(format!("budget-{budget}-gaps.txt")), &gaps).unwrap();
        std::fs::write(output.join(format!("budget-{budget}-application.txt")), format!("{application:?}\n")).unwrap();
        if let Some(baseline) = args.get(2) {
            let baseline = Path::new(baseline);
            assert_eq!(std::fs::read_to_string(baseline.join(format!("budget-{budget}-gaps.txt"))).unwrap(), gaps);
            assert!(equivalent_generated_programs(
                &std::fs::read(baseline.join(format!("budget-{budget}.bin"))).unwrap(),
                &bytes,
                Default::default(),
            ).unwrap());
        }
        println!("budget={budget} rules={} gaps={} application={:?} baseline_equivalent={}", decoded.rules().len(), found.unresolved.len(), application.status, args.len() > 2);
    }
}
