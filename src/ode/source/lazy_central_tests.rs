//! Exact behavioral equivalence against the frozen eager 977 implementation.
//! These fixtures exercise native directed ball tables and independent expected
//! bounds/errors. They do not inspect the lazy branch or count implementation calls.
use super::*;

// Frozen baseline defect_bounds, with only its receiver name made explicit.
fn eager_reference(
    source: &ExactSourceResidual,
    p: Precision,
    step: &C,
    budget: Option<(&[C], &Float)>,
) -> Result<Vec<Float>> {
    if !p.finite(step)
        || budget.is_some_and(|(values, tolerance)| {
            values.len() != source.residuals.len()
                || values.iter().any(|a| !p.finite(a))
                || !tolerance.is_finite()
                || *tolerance <= p.real(0)
        })
    {
        return Err(Error::InvalidInput(
            "exact-source residual step or budget".into(),
        ));
    }
    // A directed component one-norm is a certified upper complex radius;
    // no uncertified square root or transcendental ball operation is used.
    let radius = step
        .re
        .norm()
        .add_round(&step.im.norm(), p.bits, RoundingDirection::Up);
    let denominator_lowers = source
        .denominators
        .iter()
        .map(|denominator| denominator.lower(&radius, p))
        .collect::<Vec<_>>();
    let mut bounds = Vec::with_capacity(source.residuals.len());
    for (i, coefficients) in source.bounds.iter().enumerate() {
        let lower = &denominator_lowers[i % source.denominators.len()];
        if !lower.is_finite() || *lower <= p.real(0) {
            return Err(Error::Accuracy(
                "exact-source denominator enclosure is inconclusive on the trial disk".into(),
            ));
        }
        let mut upper = p.real(0);
        let mut width = p.real(0);
        let mut central = p.real(0);
        // All coefficients and the radius are nonnegative. Horner's rule
        // with upward rounding encloses the same polynomial disk bound.
        for (coefficient_upper, coefficient_width, coefficient_central) in coefficients.iter().rev()
        {
            upper = upper
                .mul_round(&radius, p.bits, RoundingDirection::Up)
                .add_round(coefficient_upper, p.bits, RoundingDirection::Up);
            central = central
                .mul_round(&radius, p.bits, RoundingDirection::Up)
                .add_round(coefficient_central, p.bits, RoundingDirection::Up);
            width = width
                .mul_round(&radius, p.bits, RoundingDirection::Up)
                .add_round(coefficient_width, p.bits, RoundingDirection::Up);
        }
        let bound = upper
            .div_round(lower, p.bits, RoundingDirection::Up)
            .mul_round(&radius, p.bits, RoundingDirection::Up);
        let width = width
            .div_round(lower, p.bits, RoundingDirection::Up)
            .mul_round(&radius, p.bits, RoundingDirection::Up);
        let central = central
            .div_round(lower, p.bits, RoundingDirection::Up)
            .mul_round(&radius, p.bits, RoundingDirection::Up);
        if !bound.is_finite() || !width.is_finite() {
            return Err(Error::InsufficientPrecision {
                minimum_bits: p.bits.saturating_mul(2),
                context: "nonfinite exact-source residual enclosure".into(),
            });
        }
        if let Some((values, tolerance)) = budget {
            // Use a lower component modulus for the mixed-scale budget so
            // rounding cannot enlarge the allowed local error.
            let re = values[i].re.norm();
            let im = values[i].im.norm();
            let magnitude = if re > im { re } else { im };
            let scale = if magnitude > p.real(1) {
                magnitude
            } else {
                p.real(1)
            };
            let threshold = tolerance.mul_round(&scale, p.bits, RoundingDirection::Down);
            if threshold.is_zero() {
                return Err(Error::InsufficientPrecision {
                    minimum_bits: p.bits.saturating_mul(2),
                    context: "exact-source local error budget underflows".into(),
                });
            }
            if width > threshold && width >= central {
                let ratio = width.to_rational() / threshold.to_rational();
                let extra = ratio
                    .numerator_ref()
                    .significant_bits()
                    .saturating_sub(ratio.denominator_ref().significant_bits())
                    .saturating_add(1);
                let minimum_bits = u32::try_from(extra)
                    .ok()
                    .and_then(|n| p.bits.checked_add(n))
                    .and_then(|n| n.checked_add(32))
                    .ok_or_else(|| Error::Limit("exact-source precision hint overflow".into()))?;
                return Err(Error::InsufficientPrecision {
                    minimum_bits,
                    context: format!(
                        "exact-source residual enclosure in component {i} exceeds the local error budget"
                    ),
                });
            }
        }
        bounds.push(bound);
    }
    Ok(bounds)
}

fn polynomial(p: Precision, coefficients: Vec<ComplexBall>) -> BallPolynomial {
    BallPolynomial::from_coefficients(
        &FloatField::from_rep(dyadic_ball(p, &p.zero())),
        coefficients,
        Arc::new(PolyVariable::Temporary(0)),
    )
}

fn constant_chart(p: Precision, center: i64, radius: Rational) -> ExactSourceResidual {
    let value = ComplexBall::new(
        RealBall::from_rational_ball(&Rational::from(center), &radius, p.bits),
        RealBall::exact(p.real(0)),
    );
    from_coefficients(p, vec![value])
}

fn from_coefficients(p: Precision, coefficients: Vec<ComplexBall>) -> ExactSourceResidual {
    ExactSourceResidual::from_polynomials(
        p,
        vec![polynomial(p, coefficients)],
        vec![DenominatorEnclosure {
            expanded: polynomial(p, vec![dyadic_ball(p, &p.i(1))]),
            factors: vec![],
        }],
    )
}

fn equivalent(
    chart: &ExactSourceResidual,
    p: Precision,
    step: &C,
    budget: Option<(&[C], &Float)>,
) -> Result<Vec<Float>> {
    let expected = eager_reference(chart, p, step, budget);
    let actual = chart.defect_bounds(p, step, budget);
    match (&actual, &expected) {
        (Ok(actual), Ok(expected)) => {
            assert_eq!(actual, expected);
            assert_eq!(
                actual.iter().map(Float::prec).collect::<Vec<_>>(),
                expected.iter().map(Float::prec).collect::<Vec<_>>()
            );
        }
        (Err(actual), Err(expected)) => {
            assert_eq!(
                std::mem::discriminant(actual),
                std::mem::discriminant(expected)
            );
            // Display includes every diagnostic context and precision hint.
            assert_eq!(actual.to_string(), expected.to_string());
        }
        _ => panic!("lazy/eager result mismatch: actual={actual:?}; expected={expected:?}"),
    }
    actual
}

#[test]
fn lazy_central_budget_free_complex_polynomial_matches_eager() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let coefficients = (0..17)
            .map(|k| {
                ComplexBall::from_rational_ball(
                    &Gaussian::new(Rational::from((k - 4, 7)), Rational::from((3 - k, 11))),
                    &Rational::from((1, 4096)),
                    bits,
                )
            })
            .collect();
        let chart = from_coefficients(p, coefficients);
        for step in [
            p.rational(&Rational::from((1, 3))),
            p.parse("0.125", "-0.0625").unwrap(),
        ] {
            let bounds = equivalent(&chart, p, &step, None).unwrap();
            assert!(bounds[0].is_finite() && bounds[0] > p.real(0));
        }
    }
}

#[test]
fn lazy_central_small_width_matches_eager_at_113_and_3332_bits() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let chart = constant_chart(p, 16, Rational::from((1, 1024)));
        let values = [p.i(1)];
        let tolerance = p.rational(&Rational::from((1, 16))).re;
        let step = p.rational(&Rational::from((1, 2)));
        let actual = equivalent(&chart, p, &step, Some((&values, &tolerance))).unwrap();
        assert_eq!(actual, vec![p.rational(&Rational::from((16385, 2048))).re]);
    }
}

#[test]
fn lazy_central_large_width_with_central_dominance_preserves_bound() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let chart = constant_chart(p, 16, Rational::from(4));
        let values = [p.i(1)];
        let tolerance = p.rational(&Rational::from((1, 16))).re;
        let step = p.rational(&Rational::from((1, 2)));
        let actual = equivalent(&chart, p, &step, Some((&values, &tolerance))).unwrap();
        // Width=2 exceeds the budget, but central=8 dominates it.
        assert_eq!(actual, vec![p.real(10)]);
    }
}

#[test]
fn lazy_central_width_dominance_preserves_precision_hint() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let chart = constant_chart(p, 0, Rational::from(4));
        let values = [p.i(1)];
        let tolerance = p.rational(&Rational::from((1, 16))).re;
        let step = p.rational(&Rational::from((1, 2)));
        match equivalent(&chart, p, &step, Some((&values, &tolerance))) {
            Err(Error::InsufficientPrecision {
                minimum_bits,
                context,
            }) => {
                assert_eq!(minimum_bits, bits + 38);
                assert_eq!(
                    context,
                    "exact-source residual enclosure in component 0 exceeds the local error budget"
                );
            }
            other => panic!("expected the original precision hint, got {other:?}"),
        }
    }
}

#[test]
fn lazy_central_zero_and_tiny_positive_budgets_match_eager() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let chart = constant_chart(p, 1, Rational::zero());
        let values = [p.i(1)];
        assert!(matches!(
            equivalent(&chart, p, &p.i(1), Some((&values, &p.real(0)))),
            Err(Error::InvalidInput(_))
        ));
        // Valid tolerance is positive and scale is >=1. Do not mutate global
        // MPFR exponent bounds to manufacture an unreachable zero threshold.
        let tiny = p.tolerance(10_000);
        assert!(tiny > p.real(0));
        assert_eq!(
            equivalent(&chart, p, &p.i(1), Some((&values, &tiny))).unwrap(),
            vec![p.real(1)]
        );
    }
}

#[test]
fn lazy_central_nonfinite_enclosures_match_eager() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let infinity = Float::parse("inf", Some(bits)).unwrap();
        let values = [p.i(1)];
        let tolerance = p.real(1);
        for value in [
            ComplexBall::new(
                RealBall::exact(infinity.clone()),
                RealBall::exact(p.real(0)),
            ),
            ComplexBall::new(
                RealBall::new(p.real(0), infinity.clone()),
                RealBall::exact(p.real(0)),
            ),
        ] {
            let chart = from_coefficients(p, vec![value]);
            for budget in [None, Some((values.as_slice(), &tolerance))] {
                match equivalent(&chart, p, &p.i(1), budget) {
                    Err(Error::InsufficientPrecision {
                        minimum_bits,
                        context,
                    }) => {
                        assert_eq!(minimum_bits, 2 * bits);
                        assert_eq!(context, "nonfinite exact-source residual enclosure");
                    }
                    other => panic!("expected nonfinite enclosure rejection, got {other:?}"),
                }
            }
        }
    }
}

#[test]
fn lazy_central_factored_and_inconclusive_denominators_match_eager() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let atom =
            |text| Atom::parse(text, "lazy_central_equivalence", Default::default()).unwrap();
        let variable = match atom("x").as_view() {
            AtomView::Var(value) => value.get_symbol(),
            _ => unreachable!(),
        };
        let rows = crate::ode::compile_rows(
            variable,
            &[vec![atom("1/(1-x)^10")]],
            p,
            &Default::default(),
        )
        .unwrap();
        let chart =
            ExactSourceResidual::new(p, &rows.exact_source_rows, &p.zero(), &[vec![p.i(1)]], 1)
                .unwrap();
        let step = p.rational(&Rational::from((1, 3)));
        assert_eq!(
            polynomial_lower(&chart.denominators[0].expanded, &step.re, p),
            p.real(0)
        );
        assert!(chart.denominators[0].lower(&step.re, p) > p.real(0));
        assert!(equivalent(&chart, p, &step, None).is_ok());
        let values = [p.i(1)];
        let tolerance = p.real(1000);
        assert!(equivalent(&chart, p, &step, Some((&values, &tolerance))).is_ok());
        assert!(matches!(
            equivalent(&chart, p, &p.i(2), None),
            Err(Error::Accuracy(_))
        ));
    }
}

#[test]
fn lazy_central_input_validation_matches_eager() {
    for bits in [113, 3332] {
        let p = Precision { bits };
        let chart = constant_chart(p, 1, Rational::zero());
        let infinity = p.parse("inf", "0").unwrap();
        assert!(matches!(
            equivalent(&chart, p, &infinity, None),
            Err(Error::InvalidInput(_))
        ));
        let tolerance = p.real(1);
        assert!(matches!(
            equivalent(&chart, p, &p.i(1), Some((&[], &tolerance))),
            Err(Error::InvalidInput(_))
        ));
    }
}
