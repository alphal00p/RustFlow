use crate::diffexp::{EpsilonBoundary, EpsilonSystem, transport_epsilon};
use crate::ode::SeriesSystem;
use crate::*;
use symbolica::prelude::*;

fn atom(s: &str) -> Atom {
    Atom::parse(s, "exact_source_stage", Default::default()).unwrap()
}
fn variable() -> Symbol {
    match atom("x").as_view() {
        AtomView::Var(v) => v.get_symbol(),
        _ => unreachable!(),
    }
}
fn large_system(offset: Rational) -> EpsilonSystem {
    EpsilonSystem {
        variable: variable(),
        matrices: vec![vec![
            vec![
                Atom::zero(),
                Atom::num(offset - Rational::from(Integer::from(2).pow(400))),
            ],
            vec![Atom::zero(), Atom::zero()],
        ]],
    }
}
fn boundary(p: Precision) -> EpsilonBoundary {
    EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.rational(&Integer::from(2).pow(400).into()), p.i(1)]],
    }
}
fn options(working: u32) -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: working - 20,
        series_order: 80,
        max_steps: 100,
        ..Default::default()
    }
}

#[test]
fn exact_source_normalized_cancellation_requests_precision_before_transport_admission() {
    for working in [60, 80, 100] {
        let p = Precision::decimal(working).unwrap();
        let system = large_system(Rational::one())
            .compile(p, &Default::default())
            .unwrap();
        let b = boundary(p);
        let (coefficients, chart) = system
            .local_chart(&p.zero(), &b.coefficients[0], 80, &())
            .unwrap();
        assert_eq!(coefficients[1][0], p.neg(&b.coefficients[0][0]));
        let bound = system
            .whole_segment_residual(&chart, &p.i(1))
            .unwrap()
            .unwrap();
        assert!(bound[0] >= p.real(1));
        let error = system
            .whole_segment_residual_with_budget(
                &chart,
                &p.i(1),
                &[p.zero(), p.i(1)],
                &p.tolerance(20),
            )
            .unwrap_err();
        assert!(
            matches!(error,Error::InsufficientPrecision{minimum_bits,..} if minimum_bits>p.bits)
        );
        let error = system
            .transport(
                &b,
                &[p.i(1)],
                &options(working),
                &RunContext::default(),
                true,
            )
            .unwrap_err();
        assert!(matches!(error, Error::InsufficientPrecision { .. }));
    }
    let p = Precision::decimal(120).unwrap();
    let compiled = large_system(Rational::one())
        .compile(p, &Default::default())
        .unwrap();
    let b = boundary(p);
    let (coefficients, chart) = compiled
        .local_chart(&p.zero(), &b.coefficients[0], 80, &())
        .unwrap();
    assert_eq!(
        crate::ode::evaluate_taylor(p, &coefficients, &p.i(1)).0,
        [p.i(1), p.i(1)]
    );
    assert!(
        compiled
            .whole_segment_residual(&chart, &p.i(1))
            .unwrap()
            .unwrap()
            .iter()
            .all(|bound| *bound == p.real(0))
    );
    // An exact residual does not bypass the separate endpoint perturbation
    // diagnostic: at this precision cancellation still fails that check.
    assert!(matches!(
        compiled.transport(&b, &[p.i(1)], &options(120), &RunContext::default(), true),
        Err(Error::InsufficientPrecision { .. })
    ));
    let p = Precision::decimal(160).unwrap();
    let result = large_system(Rational::one())
        .compile(p, &Default::default())
        .unwrap()
        .transport(
            &boundary(p),
            &[p.i(1)],
            &options(160),
            &RunContext::default(),
            true,
        )
        .unwrap();
    assert_eq!(result.coefficients[0], [p.i(1), p.i(1)]);
    assert_eq!(result.diagnostics.steps, 1);
}

#[test]
fn ordinary_transport_uses_the_same_exact_source_evidence() {
    let p = Precision::decimal(60).unwrap();
    let source = large_system(Rational::one());
    let compiled = DifferentialSystem {
        variable: source.variable,
        matrix: source.matrices[0].clone(),
    }
    .compile(p, &Default::default())
    .unwrap();
    let b = boundary(p);
    let error = compiled
        .transport(
            &BoundaryData {
                point: b.point,
                values: b.coefficients[0].clone(),
            },
            &[p.i(1)],
            &options(60),
            &RunContext::default(),
        )
        .unwrap_err();
    assert!(matches!(error, Error::InsufficientPrecision { .. }));
}

#[test]
fn exact_polynomial_cancellation_and_epsilon_channels_are_preserved() {
    let p = Precision::decimal(60).unwrap();
    let ordinary = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![atom("2/x")]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let (_, chart) = ordinary.local_chart(&p.i(1), &[p.i(1)], 8, &()).unwrap();
    assert_eq!(
        ordinary
            .whole_segment_residual(&chart, &p.scale(&p.i(1), 1, 2))
            .unwrap()
            .unwrap(),
        [p.real(0)]
    );
    let epsilon = EpsilonSystem {
        variable: variable(),
        matrices: vec![
            vec![vec![Atom::zero(); 2]; 2],
            vec![
                vec![Atom::zero(), Atom::num(2)],
                vec![Atom::zero(), Atom::zero()],
            ],
        ],
    }
    .compile(p, &Default::default())
    .unwrap();
    let (_, chart) = epsilon
        .local_chart(&p.zero(), &[p.zero(), p.i(3), p.zero(), p.i(7)], 8, &())
        .unwrap();
    assert_eq!(
        epsilon
            .whole_segment_residual(&chart, &p.i(1))
            .unwrap()
            .unwrap(),
        vec![p.real(0); 4]
    );
}

#[test]
fn full_sparse_source_degree_does_not_become_an_arithmetic_failure() {
    let p = Precision::decimal(60).unwrap();
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![atom("x^1000")]],
    }
    .compile(p, &Default::default())
    .unwrap();
    assert_eq!(system.exact_source_rows[0].entries[0].1.len(), 1001);
    let (_, chart) = system.local_chart(&p.zero(), &[p.i(1)], 16, &()).unwrap();
    let bound = system
        .whole_segment_residual_with_budget(&chart, &p.i(1), &[p.i(1)], &p.tolerance(50))
        .unwrap()
        .unwrap();
    assert_eq!(bound, [p.real(1)]); // True defect, zero enclosure width: ordinary step rejection.
}

#[test]
fn genuine_truncation_subdivides_and_preserves_an_analytic_result() {
    let p = Precision::decimal(45).unwrap();
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![atom("x^20")]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let result = system
        .transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.i(1)],
            &FlowOptions {
                digits: 15,
                guard_digits: 10,
                series_order: 16,
                max_steps: 3000,
                ..Default::default()
            },
            &RunContext::default(),
        )
        .unwrap();
    assert!(result.diagnostics.rejected_steps > 0);
    assert!(result.diagnostics.steps > 1);
    assert!(p.close(
        &result.values[0],
        &p.exp(&p.rational(&Rational::from((1, 21)))),
        15
    ));
}

#[test]
fn exact_complex_constants_and_dyadic_parameters_use_native_specialization() {
    let p = Precision::decimal(60).unwrap();
    let m = atom("m");
    let value = p.parse("1.125", "-0.375").unwrap();
    let exact = Complex::new(Rational::from((3, 7)), Rational::from((2, 9)));
    let expression = Atom::num(exact.clone()) - &m;
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![expression]],
    }
    .compile(p, &ahash::HashMap::from_iter([(m, value.clone())]))
    .unwrap();
    let row = &system.exact_source_rows[0];
    let quotient = row.entries[0].1[0].clone() / row.denominator[0].clone();
    assert_eq!(
        quotient,
        exact - Complex::new(value.re.to_rational(), value.im.to_rational())
    );
}

#[test]
fn lower_precision_stored_dyadics_do_not_limit_residual_arithmetic_precision() {
    let p = Precision::decimal(70).unwrap();
    let stored = Precision::decimal(5).unwrap().parse("0.1", "0").unwrap();
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![Atom::one()]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let (coefficients, chart) = system
        .local_chart(&p.zero(), std::slice::from_ref(&stored), 80, &())
        .unwrap();
    let (values, _) = crate::ode::evaluate_taylor(p, &coefficients, &p.i(1));
    let bounds = system
        .whole_segment_residual_with_budget(&chart, &p.i(1), &values, &p.tolerance(60))
        .unwrap()
        .unwrap();
    assert!(bounds[0] <= p.tolerance(60));
    // This checks the residual of the polynomial with exact stored dyadic
    // coefficients. Endpoint conditioning may separately reject the low-bit
    // boundary; this test must not bypass or assert that policy away.
    let exact_dyadic = p.rational(&stored.re.to_rational());
    assert!(p.close(&values[0], &p.mul(&exact_dyadic, &p.exp(&p.i(1))), 20));
}

#[test]
fn rational_mobius_pullback_retains_exact_mapped_rows() {
    let p = Precision::decimal(70).unwrap();
    let z = variable();
    let x = atom("u");
    let map = (Atom::one() + Atom::var(z)) / (Atom::one() - Atom::var(z));
    let source = Atom::one() / &x;
    let pulled = (source.replace(x.clone()).with(map.clone()) * map.derivative(z))
        .together()
        .cancel();
    let system = DifferentialSystem {
        variable: z,
        matrix: vec![vec![pulled]],
    }
    .compile(p, &Default::default())
    .unwrap();
    assert_eq!(system.exact_source_rows[0].denominator.len(), 3);
    assert_eq!(system.exact_source_rows[0].entries[0].1.len(), 1);
    let quarter = p.rational(&Rational::from((1, 4)));
    let result = system
        .transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[quarter],
            &options(70),
            &RunContext::default(),
        )
        .unwrap();
    assert!(p.close(&result.values[0], &p.rational(&Rational::from((5, 3))), 20));
}

#[test]
fn denominator_enclosure_is_inconclusive_instead_of_requesting_precision() {
    let p = Precision::decimal(60).unwrap();
    let system = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![atom("1/(1-x)")]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let (_, chart) = system.local_chart(&p.zero(), &[p.i(1)], 8, &()).unwrap();
    let error = system
        .whole_segment_residual_with_budget(&chart, &p.i(1), &[p.i(1)], &p.tolerance(20))
        .unwrap_err();
    assert!(matches!(error, Error::Accuracy(_)));
}

#[test]
fn nondyadic_severe_cancellation_fails_boundedly_without_false_metadata() {
    // Under the existing working-derived local tolerance, source roundoff and
    // tolerance tighten together. Preserve that conservative limitation.
    for working in [120, 140, 160] {
        let p = Precision::decimal(working).unwrap();
        let error = large_system(Rational::from((1, 3)))
            .compile(p, &Default::default())
            .unwrap()
            .transport(
                &boundary(p),
                &[p.i(1)],
                &options(working),
                &RunContext::default(),
                true,
            )
            .unwrap_err();
        assert!(matches!(error, Error::InsufficientPrecision { .. }));
    }
    // The shared refinement wrapper retries this arithmetic failure, but the
    // moving local tolerance remains unresolved within its finite budget.
    let error = transport_epsilon(
        &large_system(Rational::from((1, 3))),
        |p| Ok(boundary(p)),
        &[atom("1")],
        &options(120),
        &RunContext::default(),
        true,
    )
    .unwrap_err();
    assert!(matches!(error, Error::InsufficientPrecision { .. }));
}
