use crate::{
    diffexp::{EpsilonBoundary, EpsilonSystem, transport_epsilon},
    ode::{SeriesSystem, compile_rows},
    *,
};
use symbolica::prelude::*;
fn atom(s: &str) -> Atom {
    Atom::parse(s, "whole_residual_tests", Default::default()).unwrap()
}
fn variable() -> Symbol {
    match atom("x").as_view() {
        AtomView::Var(v) => v.get_symbol(),
        _ => unreachable!(),
    }
}
fn hidden() -> EpsilonSystem {
    EpsilonSystem {
        variable: variable(),
        matrices: vec![
            vec![vec![atom("0")]],
            vec![vec![atom("2*1001*1002*1003*x^1000*(x-1/2)*(x-1)")]],
        ],
    }
}
#[test]
fn sparse_alias_full_step_rejected_with_exact_zero_point_defects() {
    let p = Precision::decimal(60).unwrap();
    let sys = hidden().compile(p, &Default::default()).unwrap();
    let (coefficients, chart) = sys
        .local_chart(&p.zero(), &[p.i(1), p.zero()], 80, &())
        .unwrap();
    assert!(
        coefficients
            .iter()
            .skip(1)
            .flatten()
            .all(|v| *v == p.zero())
    );
    let bounds = sys
        .whole_segment_residual(&chart, &p.i(1))
        .unwrap()
        .unwrap();
    assert_eq!(bounds[0], p.real(0));
    assert!(bounds[1] > p.real(999));
}
#[test]
fn hidden_system_never_receives_false_verified_digits() {
    for strategy in [StepSizeStrategy::Halving, StepSizeStrategy::Bracketed] {
        let options = FlowOptions {
            digits: 20,
            guard_digits: 40,
            series_order: 80,
            max_steps: 1000,
            step_size_strategy: strategy,
            ..Default::default()
        };
        let answer = transport_epsilon(
            &hidden(),
            |p| {
                Ok(EpsilonBoundary {
                    point: p.zero(),
                    leading: 0,
                    coefficients: vec![vec![p.i(1)], vec![p.zero()]],
                })
            },
            &[atom("1")],
            &options,
            &RunContext::default(),
            true,
        );
        match answer {
            Ok(solution) => {
                let p = Precision::decimal(80).unwrap();
                assert_eq!(solution.verified_digits, Some(20));
                assert!(p.close(&solution.coefficients[1][0], &p.i(-999), 20));
                assert!(solution.diagnostics.steps > 1);
                eprintln!(
                    "hidden {strategy:?}: {:?}, epsilon1={}",
                    solution.diagnostics, solution.coefficients[1][0]
                );
            }
            Err(Error::Accuracy(_) | Error::Limit(_)) => {}
            Err(error) => panic!("unexpected error {error}"),
        }
    }
}
#[test]
fn exact_matrix_polynomial_and_all_channels_cancel() {
    // Y0=[1,2], Y1=[x,x], Y2=[x^2/4,x^2/2] for A1=[[0,1/2],[1,0]].
    let p = Precision::decimal(60).unwrap();
    let x = variable();
    let rows = compile_rows(
        x,
        &[
            vec![atom("0"), atom("0"), atom("0"), atom("1/2")],
            vec![atom("0"), atom("0"), atom("1"), atom("0")],
        ],
        p,
        &Default::default(),
    )
    .unwrap();
    let coefficients = vec![
        vec![p.i(1), p.i(2), p.zero(), p.zero(), p.zero(), p.zero()],
        vec![p.zero(), p.zero(), p.i(1), p.i(1), p.zero(), p.zero()],
        vec![
            p.zero(),
            p.zero(),
            p.zero(),
            p.zero(),
            p.rational(&Rational::from((1, 4))),
            p.rational(&Rational::from((1, 2))),
        ],
    ];
    let chart = ode::residual::RationalResidualChart::new(
        p,
        &rows.polynomial_rows,
        &p.zero(),
        &coefficients,
        3,
    )
    .unwrap();
    assert!(
        chart
            .defect_bounds(p, &p.i(3))
            .unwrap()
            .iter()
            .all(|v| *v == p.real(0))
    );
}
#[test]
fn rational_denominator_cancellation_and_inconclusive_disk() {
    // y=(1-x)^-1, P=1+x+...+x^8 has residual -9*x^8/(1-x).
    let p = Precision::decimal(60).unwrap();
    let rows = compile_rows(variable(), &[vec![atom("1/(1-x)")]], p, &Default::default()).unwrap();
    let coefficients = vec![vec![p.i(1)]; 9];
    let chart = ode::residual::RationalResidualChart::new(
        p,
        &rows.polynomial_rows,
        &p.zero(),
        &coefficients,
        1,
    )
    .unwrap();
    let h = p.rational(&Rational::from((1, 4)));
    let bound = chart.defect_bounds(p, &h).unwrap();
    let expected = p.rational(&Rational::from((12, 4_i64.pow(9))));
    assert!(p.close(
        &ComplexFloat::new(bound[0].clone(), p.real(0)),
        &expected,
        50
    ));
    assert!(matches!(
        chart.defect_bounds(p, &p.i(1)),
        Err(Error::Accuracy(_))
    ));
}
#[test]
fn complex_center_and_step_native_shift() {
    let p = Precision::decimal(60).unwrap();
    let rows = compile_rows(variable(), &[vec![atom("2*x")]], p, &Default::default()).unwrap();
    let center = p.parse("0", "1").unwrap();
    let coefficients = vec![vec![p.i(1)]];
    let chart = ode::residual::RationalResidualChart::new(
        p,
        &rows.polynomial_rows,
        &center,
        &coefficients,
        1,
    )
    .unwrap();
    let step = p.parse("0", "0.25").unwrap();
    let bound = chart.defect_bounds(p, &step).unwrap();
    assert_eq!(bound[0], p.rational(&Rational::from((5, 8))).re);
}
