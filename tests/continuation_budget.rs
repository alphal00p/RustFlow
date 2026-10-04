use symbolica::prelude::*;
use symbolica_amflow::{
    BoundaryData, DifferentialSystem, Error, FlowOptions, Precision, Result, RunContext,
};

fn attempt(digits: u32, guard_digits: u32, max_steps: usize) -> Result<()> {
    let p = Precision::decimal(digits + guard_digits)?;
    let system = DifferentialSystem {
        variable: symbol!("continuation_budget::x"),
        matrix: vec![vec![Atom::one()]],
    };
    system.compile(p, &Default::default())?.transport(
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.i(1)],
        &FlowOptions {
            digits,
            guard_digits,
            series_order: 8,
            max_steps,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    Ok(())
}

#[test]
fn rejected_step_exhaustion_is_a_limit_not_an_accuracy_failure() {
    match attempt(20, 20, 1) {
        Err(Error::Limit(message)) => assert!(message.contains("0 accepted, 1 rejected")),
        other => panic!("expected exhausted continuation budget, got {other:?}"),
    }
}

#[test]
fn unresolved_tail_with_available_budget_remains_an_accuracy_failure() {
    match attempt(300, 50, 1000) {
        Err(Error::Accuracy(message)) => assert!(message.contains("32 local rejections")),
        other => panic!("expected a truncation accuracy failure, got {other:?}"),
    }
}
