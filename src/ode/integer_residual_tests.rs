use super::*;
use crate::ode::{DifferentialSystem, compile_rows};
use symbolica::domains::float::RealBall;

fn atom(a: &str) -> Atom {
    Atom::parse(a, "adaptive_integer_tests", Default::default()).unwrap()
}
fn variable() -> Symbol {
    match atom("x").as_view() {
        AtomView::Var(a) => a.get_symbol(),
        _ => unreachable!(),
    }
}
fn rows(p: Precision, matrix: &[&[&str]]) -> Vec<ExactPolynomialRow> {
    compile_rows(
        variable(),
        &matrix
            .iter()
            .map(|r| r.iter().map(|a| atom(a)).collect())
            .collect::<Vec<_>>(),
        p,
        &Default::default(),
    )
    .unwrap()
    .exact_source_rows
}
fn integer(a: &AdaptiveResidual) -> (&ExactSourceResidual, &Arc<LazyBall>) {
    match a {
        AdaptiveResidual::Integer { source, fallback } => (source, fallback),
        _ => panic!("expected native integer selection"),
    }
}
fn contains(a: &RealBall, q: &Rational) -> bool {
    a.lower_bound().to_rational() <= *q && *q <= a.upper_bound().to_rational()
}
// Independent native Gaussian-rational reference; no shared-denominator
// arithmetic or production IntegerPair helper is used here.
fn reference(
    rows: &[ExactPolynomialRow],
    center: &C,
    coefficients: &[Vec<C>],
    channels: usize,
) -> Vec<ExactPolynomial> {
    let n = rows.len();
    let ring = FloatField::from_rep(Gaussian::new(Rational::zero(), Rational::zero()));
    let polynomial =
        |a| ExactPolynomial::from_coefficients(&ring, a, Arc::new(PolyVariable::Temporary(0)));
    let y = (0..n * channels)
        .map(|i| polynomial(coefficients.iter().map(|r| gaussian(&r[i])).collect()))
        .collect::<Vec<_>>();
    let center = gaussian(center);
    let mut result = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let d = polynomial(row.denominator.clone()).shift_var(&center);
        for channel in 0..channels {
            let mut residual = &d * &y[channel * n + i].derivative();
            for (column, a) in &row.entries {
                let shift = column / n;
                if shift <= channel {
                    residual = &residual
                        - &(&polynomial(a.clone()).shift_var(&center)
                            * &y[(channel - shift) * n + column % n]);
                }
            }
            result.push(residual);
        }
        result.push(d);
        result.extend(
            row.denominator_factors
                .iter()
                .map(|(a, _)| polynomial(a.clone()).shift_var(&center)),
        );
    }
    result
}

#[test]
fn integer_coefficients_equal_gaussian_reference_for_complex_epsilon_channels() {
    for bits in [64, 216, 548] {
        let p = Precision { bits };
        let rows = rows(
            p,
            &[
                &["(1+x)/(3-2*x)", "2/(5+x)", "0", "x/(7-3*x)"],
                &["1/(2+x)", "-1/(3-2*x)", "3/(5+x)", "0"],
            ],
        );
        // Deliberately distinct binary scales across six channels/components;
        // source rows have distinct odd factors. The stored precision exceeds
        // the arithmetic precision and must not be rounded before conversion.
        let storage = Precision { bits: bits + 97 };
        let coefficients = (0..5)
            .map(|k| {
                (0..6)
                    .map(|i| {
                        storage
                            .parse(
                                &format!("{}.{}", i + 1, 173 + k * 11 + i),
                                &format!("-0.{}", 197 + k * 7 + i),
                            )
                            .unwrap()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let center = storage.parse("0.137", "-0.219").unwrap();
        let expected = reference(&rows, &center, &coefficients, 3);
        let cancellation = CancellationToken::default();
        let mut budget = Budget {
            limits: Limits::default(),
            work: 0,
            bytes: 0,
            cancellation: &cancellation,
        };
        let mut exact = Vec::new();
        let new = construct_observed(p, &rows, &center, &coefficients, 3, &mut budget, |a| {
            exact.push(a.clone())
        })
        .unwrap();
        assert_eq!(
            exact, expected,
            "stored-dyadic exact equality at {bits} bits"
        );
        let old = ExactSourceResidual::new(p, &rows, &center, &coefficients, 3).unwrap();
        let mut index = 0;
        for i in 0..rows.len() {
            let mut compare = |a: &BallPolynomial, b: &BallPolynomial| {
                let q = &expected[index];
                index += 1;
                let zero = Gaussian::new(Rational::zero(), Rational::zero());
                let zero_ball = exact_ball(&zero, p);
                for k in 0..a
                    .coefficients()
                    .len()
                    .max(b.coefficients().len())
                    .max(q.coefficients().len())
                {
                    let a = a.coefficients().get(k).unwrap_or(&zero_ball);
                    let b = b.coefficients().get(k).unwrap_or(&zero_ball);
                    let q = q.coefficients().get(k).unwrap_or(&zero);
                    assert!(contains(&a.re, &q.re) && contains(&b.re, &q.re));
                    assert!(contains(&a.im, &q.im) && contains(&b.im, &q.im));
                }
            };
            for channel in 0..3 {
                compare(
                    &old.residuals[channel * rows.len() + i],
                    &new.residuals[channel * rows.len() + i],
                );
            }
            compare(&old.denominators[i].expanded, &new.denominators[i].expanded);
            for ((a, _), (b, _)) in old.denominators[i]
                .factors
                .iter()
                .zip(&new.denominators[i].factors)
            {
                compare(a, b);
            }
        }
        assert_eq!(index, expected.len());
    }
}

#[test]
fn tighter_passing_defect_does_not_construct_old_balls() {
    let p = Precision { bits: 216 };
    let rows = rows(p, &[&["1/3"]]);
    let values = vec![vec![p.i(1)]];
    let chart =
        AdaptiveResidual::new(p, &rows, &p.zero(), &values, 1, &Default::default()).unwrap();
    let (_, fallback) = integer(&chart);
    let step = p.rational(&Rational::from((1, 1000)));
    let bound = chart
        .defect_bounds(
            p,
            &step,
            Some((&values[0], &p.rational(&Rational::from((1, 100))).re)),
        )
        .unwrap();
    assert!(bound[0] > p.real(0));
    assert!(fallback.residual.lock().unwrap().is_none());
}

#[test]
fn rejected_exact_defect_recovers_precision_hint_and_reuses_ball_owner() {
    for bits in [216, 415, 548] {
        let p = Precision { bits };
        let large = Rational::from(Integer::from(2).pow(400));
        let system = DifferentialSystem {
            variable: variable(),
            matrix: vec![
                vec![Atom::zero(), Atom::num(Rational::from((1, 3)) - &large)],
                vec![Atom::zero(), Atom::zero()],
            ],
        }
        .compile(p, &Default::default())
        .unwrap();
        let coefficients = system
            .taylor(&p.zero(), &[p.rational(&large), p.i(1)], 80)
            .unwrap();
        let chart = AdaptiveResidual::new(
            p,
            &system.exact_source_rows,
            &p.zero(),
            &coefficients,
            1,
            &Default::default(),
        )
        .unwrap();
        let (source, fallback) = integer(&chart);
        assert!(fallback.residual.lock().unwrap().is_none());
        let step = p.i(1);
        let values = crate::ode::evaluate_taylor(p, &coefficients, &step).0;
        let tolerance = p.tolerance(bits / 4);
        let bounds = source.defect_bounds(p, &step, None).unwrap();
        assert!(bounds[0] > tolerance);
        let expected =
            ExactSourceResidual::new(p, &system.exact_source_rows, &p.zero(), &coefficients, 1)
                .unwrap()
                .defect_bounds(p, &step, Some((&values, &tolerance)))
                .unwrap_err();
        assert!(matches!(expected, Error::InsufficientPrecision { .. }));
        let result = chart
            .defect_bounds(p, &step, Some((&values, &tolerance)))
            .unwrap_err();
        assert_eq!(result.to_string(), expected.to_string());
        let owner = fallback
            .residual
            .lock()
            .unwrap()
            .as_ref()
            .unwrap()
            .residuals
            .as_ptr();
        let cloned = chart.clone();
        let (_, other) = integer(&cloned);
        assert!(Arc::ptr_eq(fallback, other));
        assert!(matches!(
            cloned.defect_bounds(p, &step, Some((&values, &tolerance))),
            Err(Error::InsufficientPrecision { .. })
        ));
        assert_eq!(
            owner,
            other
                .residual
                .lock()
                .unwrap()
                .as_ref()
                .unwrap()
                .residuals
                .as_ptr()
        );
    }
}

#[test]
fn metadata_and_cumulative_limits_fall_back_without_changing_support() {
    let p = Precision::decimal(60).unwrap();
    let cancellation = CancellationToken::default();
    let values = vec![vec![p.i(1)]];
    let ordinary = rows(p, &[&["1/(1-x)^64"]]);
    assert!(matches!(
        AdaptiveResidual::new(
            p,
            &ordinary,
            &p.rational(&Rational::from((1, 10))),
            &values,
            1,
            &cancellation
        )
        .unwrap(),
        AdaptiveResidual::Ball { .. }
    ));
    let sparse = rows(p, &[&["1/(1-x^257)"]]);
    assert!(matches!(
        AdaptiveResidual::new(p, &sparse, &p.zero(), &values, 1, &cancellation).unwrap(),
        AdaptiveResidual::Integer { .. }
    ));
    let simple = rows(p, &[&["1/(1-x)"]]);
    for limits in [
        Limits {
            work: 1,
            ..Default::default()
        },
        Limits {
            bytes: 1,
            ..Default::default()
        },
        Limits {
            integer_bits: 8,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            AdaptiveResidual::with_limits(p, &simple, &p.zero(), &values, 1, &cancellation, limits)
                .unwrap(),
            AdaptiveResidual::Ball { .. }
        ));
    }
    let mut huge = simple.clone();
    huge[0].entries[0].1[0].re = Rational::from(Integer::from(2).pow(5000));
    assert!(matches!(
        AdaptiveResidual::new(p, &huge, &p.zero(), &values, 1, &cancellation).unwrap(),
        AdaptiveResidual::Ball { .. }
    ));
    // A modest-precision Float with a huge binary exponent is rejected BEFORE
    // to_rational; the existing ball constructor can represent it compactly.
    for reciprocal in [false, true] {
        let large = Float::with_val(p.bits, 2).pow(100000);
        let a = if reciprocal { p.real(1) / large } else { large };
        let values = vec![vec![C::new(a, p.real(0))]];
        assert!(matches!(
            AdaptiveResidual::new(p, &simple, &p.zero(), &values, 1, &cancellation).unwrap(),
            AdaptiveResidual::Ball { .. }
        ));
    }
}

#[test]
fn cancellation_and_invalid_inputs_never_initialize_or_poison_fallback() {
    let p = Precision::decimal(40).unwrap();
    let rows = rows(p, &[&["1/(1-x)"]]);
    let token = CancellationToken::default();
    let chart = AdaptiveResidual::new(p, &rows, &p.zero(), &[vec![p.i(1)]], 1, &token).unwrap();
    let (_, fallback) = integer(&chart);
    assert!(matches!(
        chart.defect_bounds(p, &p.i(1), Some((&[], &p.real(1)))),
        Err(Error::InvalidInput(_))
    ));
    assert!(fallback.residual.lock().unwrap().is_none());
    token.cancel();
    assert!(matches!(
        chart.defect_bounds(p, &p.i(1), None),
        Err(Error::Cancelled)
    ));
    assert!(fallback.residual.lock().unwrap().is_none());
    assert!(matches!(
        AdaptiveResidual::new(p, &rows, &p.zero(), &[vec![p.i(1)]], 1, &token),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        AdaptiveResidual::new(p, &rows, &p.zero(), &[vec![]], 1, &Default::default()),
        Err(Error::InvalidInput(_))
    ));
    // Direct owner cancellation is checked before taking its lock as well.
    let guard = fallback.residual.lock().unwrap();
    assert!(matches!(
        fallback.bounds(&p.i(1), None),
        Err(Error::Cancelled)
    ));
    drop(guard);
}

#[test]
fn factorized_denominator_holes_are_not_removed_by_exact_cancellation() {
    let p = Precision::decimal(60).unwrap();
    let rows = rows(p, &[&["1/(1-x)^10"]]);
    let values = vec![vec![p.zero()]];
    let chart =
        AdaptiveResidual::new(p, &rows, &p.zero(), &values, 1, &Default::default()).unwrap();
    let (_, fallback) = integer(&chart);
    assert_eq!(
        chart
            .defect_bounds(p, &p.rational(&Rational::from((1, 3))), None)
            .unwrap(),
        vec![p.real(0)]
    );
    assert!(matches!(
        chart.defect_bounds(p, &p.i(2), Some((&values[0], &p.real(1)))),
        Err(Error::Accuracy(_))
    ));
    assert!(fallback.residual.lock().unwrap().is_some());
}

#[test]
fn public_options_reach_rational_and_epsilon_mapped_chart_owners() {
    use crate::ode::SeriesSystem;
    use crate::{FlowOptions, ResidualArithmetic, RunContext};
    let p = Precision::decimal(60).unwrap();
    let source = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![atom("(1+eps)/(1-x)")]],
    };
    let epsilon = match atom("eps").as_view() {
        AtomView::Var(a) => a.get_symbol(),
        _ => unreachable!(),
    };
    let ordinary = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![atom("1/(1-x)")]],
    }
    .compile(p, &Default::default())
    .unwrap();
    let epsilon = crate::diffexp::EpsilonSystem::from_differential_system(&source, epsilon, 3)
        .unwrap()
        .compile(p, &Default::default())
        .unwrap();
    let center = p.zero();
    let coordinate = crate::ode::TaylorCoordinate::balanced(
        p,
        &center,
        &p.rational(&Rational::from((1, 4))),
        &ordinary.poles,
    )
    .unwrap();
    for mode in [
        ResidualArithmetic::Ball,
        ResidualArithmetic::AdaptiveInteger,
    ] {
        let options = FlowOptions {
            residual_arithmetic: mode,
            ..Default::default()
        };
        let context = RunContext::default();
        let (_, chart) = ordinary
            .local_chart_with_options(&center, &[p.i(1)], 16, &(), &options, &context)
            .unwrap();
        assert_eq!(chart.selected_arithmetic(), Some(mode));
        let (_, mapped, _) = ordinary
            .mapped_chart_with_options(&coordinate, &center, &[p.i(1)], 16, &(), &options, &context)
            .unwrap();
        assert_eq!(mapped.selected_arithmetic(), Some(mode));
        let values = [p.i(1), p.zero(), p.zero(), p.zero()];
        let (_, chart) = epsilon
            .local_chart_with_options(&center, &values, 16, &(), &options, &context)
            .unwrap();
        assert_eq!(chart.selected_arithmetic(), Some(mode));
        let (_, mapped, _) = epsilon
            .mapped_chart_with_options(&coordinate, &center, &values, 16, &(), &options, &context)
            .unwrap();
        assert_eq!(mapped.selected_arithmetic(), Some(mode));
        context.cancellation.cancel();
        assert!(matches!(
            ordinary.local_chart_with_options(&center, &[p.i(1)], 16, &(), &options, &context),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            epsilon.mapped_chart_with_options(
                &coordinate,
                &center,
                &values,
                16,
                &(),
                &options,
                &context
            ),
            Err(Error::Cancelled)
        ));
    }
    let (_, plain) = ordinary.local_chart(&center, &[p.i(1)], 16, &()).unwrap();
    assert_eq!(plain.selected_arithmetic(), Some(ResidualArithmetic::Ball));
}
