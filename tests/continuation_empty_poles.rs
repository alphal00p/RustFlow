use symbolica::prelude::*;
use symbolica_amflow::{
    BoundaryData, DifferentialSystem, FlowOptions, Precision, Result, RunContext,
};
#[test]
fn entire_system_reaches_a_real_endpoint_after_complex_waypoints() -> Result<()> {
    let p = Precision::decimal(50)?;
    let system = DifferentialSystem {
        variable: symbol!("entire_contour::x"),
        matrix: vec![vec![Atom::new()]],
    };
    let result = system.compile(p, &Default::default())?.transport(
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.complex(0, 1), p.i(1)],
        &FlowOptions {
            digits: 20,
            guard_digits: 30,
            series_order: 16,
            max_steps: 4,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    assert_eq!(result.point, p.i(1));
    assert_eq!(result.values, vec![p.i(1)]);
    assert_eq!(result.diagnostics.steps, 2);
    Ok(())
}

#[test]
fn entire_nonzero_kernel_is_checked_across_complex_waypoints() -> Result<()> {
    let p = Precision::decimal(50)?;
    let system = DifferentialSystem {
        variable: symbol!("entire_contour::x"),
        matrix: vec![vec![Atom::one()]],
    };
    let result = system.compile(p, &Default::default())?.transport(
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.complex(0, 1), p.i(1)],
        &FlowOptions {
            digits: 20,
            guard_digits: 30,
            series_order: 64,
            max_steps: 4,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    assert!(p.close(&result.values[0], &p.exp(&p.i(1)), 40));
    assert_eq!(result.point, p.i(1));
    assert_eq!(result.diagnostics.steps, 2);
    Ok(())
}

#[test]
fn entire_polynomial_preserves_accuracy_across_complex_direction_changes() -> Result<()> {
    let p = Precision::decimal(60)?;
    let x = symbol!("entire_step_proposal::x");
    let system = DifferentialSystem {
        variable: x,
        matrix: vec![vec![Atom::var(x).pow(2) + Atom::one()]],
    }
    .compile(p, &Default::default())?;
    let result = system.transport(
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.complex(0, 1), p.i(1), p.complex(0, -1)],
        &FlowOptions {
            digits: 20,
            guard_digits: 30,
            series_order: 24,
            max_steps: 1000,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    let expected = p.exp(&p.scale(&p.complex(0, -2), 1, 3));
    assert!(p.close(&result.values[0], &expected, 35));
    assert_eq!(result.point, p.complex(0, -1));
    assert!(result.diagnostics.rejected_steps > 0);
    Ok(())
}

#[test]
fn entire_exponential_reverses_direction_after_rejected_steps() -> Result<()> {
    let p = Precision::decimal(60)?;
    let system = DifferentialSystem {
        variable: symbol!("entire_step_exponential::x"),
        matrix: vec![vec![Atom::one()]],
    }
    .compile(p, &Default::default())?;
    let result = system.transport(
        &BoundaryData {
            point: p.zero(),
            values: vec![p.i(1)],
        },
        &[p.complex(1, 1), p.complex(-1, 1), p.zero()],
        &FlowOptions {
            digits: 20,
            guard_digits: 20,
            series_order: 24,
            max_steps: 1000,
            ..Default::default()
        },
        &RunContext::default(),
    )?;
    assert_eq!(result.point, p.zero());
    assert!(p.close(&result.values[0], &p.i(1), 25));
    assert!(result.diagnostics.rejected_steps > 0);
    Ok(())
}
