//! Exact-specialization pole regressions and source/enclosure invariants.
use super::*;

fn variable() -> Symbol {
    symbol!("repeated_poles_tests::t")
}

fn parameter(name: &str) -> Atom {
    Atom::var(symbol!(name))
}

/// Compare formal factors without imposing J=i or numerical parameters.
fn formal_rational(expression: &Atom) -> Result<RationalPolynomial<IntegerRing, u16>> {
    let variables = Arc::new(vec![
        PolyVariable::Symbol(crate::family::imaginary_parameter()),
        PolyVariable::Symbol(variable()),
    ]);
    expression
        .try_to_rational_polynomial(&Q, &Z, Some(variables))
        .map_err(|error| Error::InvalidInput(error.to_string()))
}

fn exact_polynomial(coefficients: Vec<source::Gaussian>) -> ExactPolePolynomial {
    ExactPolePolynomial::from_coefficients(
        &FloatField::from_rep(source::Gaussian::new(Rational::zero(), Rational::zero())),
        coefficients,
        Arc::new(PolyVariable::Symbol(variable())),
    )
}

fn specialize(
    expression: &Atom,
    p: Precision,
    values: &ahash::HashMap<Atom, C>,
) -> Result<ExactPolePolynomial> {
    Ok(exact_polynomial(
        source::ExactSpecialization::new(p, values)?.polynomial(expression, variable())?,
    ))
}

fn square_free(
    expression: &Atom,
    p: Precision,
    values: &ahash::HashMap<Atom, C>,
) -> Result<ExactPolePolynomial> {
    square_free_pole_polynomial(
        source::ExactSpecialization::new(p, values)?.polynomial(expression, variable())?,
        variable(),
    )
}

/// Expanded factors are irreducible over free J, but repeated over Q(i).
fn quadratic_cases() -> [(Atom, (i64, i64, i64), (usize, usize)); 4] {
    let t = Atom::var(variable());
    let j = Atom::var(crate::family::imaginary_parameter());
    let jm = &j - 1;
    let b = &t + Atom::num(Rational::from((5, 2))) * &jm;
    let c = &t + &jm;
    [
        (
            Atom::num(2) * t.clone().pow(2) + Atom::num(10) * &jm * &t - Atom::num(25) * &j,
            (5, -5, 2),
            (2, 1),
        ),
        (Atom::num(2) * b.pow(2), (5, -5, 2), (1, 2)),
        (
            t.pow(2) + Atom::num(2) * &jm * &t - Atom::num(2) * &j,
            (1, -1, 1),
            (2, 1),
        ),
        (c.pow(2), (1, -1, 1), (1, 2)),
    ]
}

fn quartic() -> Atom {
    let t = Atom::var(variable());
    let j = Atom::var(crate::family::imaginary_parameter());
    Atom::num(16) * t.clone().pow(4) + Atom::num(74) * (&j - 1) * t.clone().pow(3)
        - Atom::num(216) * &j * t.clone().pow(2)
        + Atom::num(125) * (Atom::one() + &j) * &t
        - 50
}

#[test]
fn public_expanded_and_factored_quadratics_compile_at_482_bits() -> Result<()> {
    let p = Precision { bits: 482 };
    let mut found = Vec::new();
    for (denominator, (re, im, divisor), _) in quadratic_cases() {
        let system = DifferentialSystem {
            variable: variable(),
            matrix: vec![vec![Atom::one() / denominator]],
        };
        let compiled = system.compile(p, &Default::default())?;
        assert_eq!(compiled.poles.len(), 1);
        assert_eq!(compiled.poles[0], p.scale(&p.complex(re, im), 1, divisor));
        assert_eq!(compiled.source.rows, system.matrix);
        // Multiplicity reduction must not replace the quadratic row denominator.
        assert_eq!(compiled.matrix[0][0].denominator.len(), 3);
        assert_eq!(compiled.exact_source_rows[0].denominator.len(), 3);
        found.push(compiled.poles[0].clone());
    }
    assert_eq!(found[0], found[1]);
    assert_eq!(found[2], found[3]);
    Ok(())
}

#[test]
fn exact_powers_and_gaussian_content_leave_only_distinct_poles() -> Result<()> {
    let p = Precision { bits: 482 };
    let t = Atom::var(variable());
    let content = Atom::num(Complex::new(Rational::from(3), Rational::from(2)));
    let original = &content * (&t - 2).pow(3) * (&t + 1).pow(2);
    let exact = specialize(&original, p, &Default::default())?;
    assert_eq!(exact.degree(), 5);
    assert_eq!(
        exact.lcoeff(),
        source::Gaussian::new(Rational::from(3), Rational::from(2))
    );
    let expected = specialize(&((&t - 2) * (&t + 1)), p, &Default::default())?;
    let actual = square_free(&original, p, &Default::default())?;
    assert_eq!(actual, expected);
    assert_eq!(actual.degree(), 2);
    assert!(actual.gcd_euclidean(&actual.derivative()).is_one());
    assert_eq!(
        exact,
        &specialize(&(&t - 2), p, &Default::default())?.pow(3)
            * &specialize(&(&t + 1), p, &Default::default())?
                .pow(2)
                .mul_coeff(&source::Gaussian::new(Rational::from(3), Rational::from(2)))
    );
    Ok(())
}

#[test]
fn later_quartic_stays_degree_four_and_square_free() -> Result<()> {
    let p = Precision { bits: 482 };
    let t = Atom::var(variable());
    let j = Atom::var(crate::family::imaginary_parameter());
    let jm = &j - 1;
    let quartic = quartic();
    let exact = specialize(&quartic, p, &Default::default())?;
    let actual = square_free(&quartic, p, &Default::default())?;
    assert_eq!(exact.degree(), 4);
    assert!(exact.gcd_euclidean(&exact.derivative()).is_one());
    assert_eq!(actual, exact.make_monic());
    assert_eq!(actual.degree(), 4);
    let product = (&t + Atom::num(Rational::from((1, 2))) * &jm)
        * (&t + Atom::num(Rational::from((5, 2))) * &jm)
        * (&t + &jm)
        * (&t + Atom::num(Rational::from((5, 8))) * &jm);
    assert_eq!(actual, specialize(&product, p, &Default::default())?);
    Ok(())
}

#[test]
fn every_compile_recomputes_repetition_and_exact_degree_collapse() -> Result<()> {
    let p = Precision { bits: 482 };
    let t = Atom::var(variable());
    let a = parameter("repeated_poles_tests::a");
    let b = parameter("repeated_poles_tests::b");
    let c = parameter("repeated_poles_tests::c");
    let factor = &a * t.pow(2) + &b * &t + &c;
    let context = RunContext::default();
    let prepared = PreparedRows::new(variable(), &[vec![Atom::one() / &factor]], &context)?;
    for (av, bv, cv, expected_root) in [
        (1, -2, 1, Some(1)),
        (0, 2, -4, Some(2)),
        (0, 0, 3, None),
        (0, 0, 0, None),
        (1, -2, 1, Some(1)),
    ] {
        let values = [
            (a.clone(), p.i(av)),
            (b.clone(), p.i(bv)),
            (c.clone(), p.i(cv)),
        ]
        .into_iter()
        .collect();
        let compiled = prepared.compile(p, &values, &context);
        if av == 0 && bv == 0 && cv == 0 {
            assert!(matches!(compiled, Err(Error::Numerical(_))));
            assert!(matches!(
                square_free(&factor, p, &values),
                Err(Error::Numerical(_))
            ));
        } else {
            let compiled = compiled?;
            let reduced = square_free(&factor, p, &values)?;
            assert_eq!(reduced.degree(), usize::from(expected_root.is_some()));
            assert_eq!(
                compiled.poles,
                expected_root.map(|r| vec![p.i(r)]).unwrap_or_default()
            );
            if expected_root.is_none() {
                assert_eq!(compiled.exact_source_rows[0].denominator.len(), 1);
                assert!(!compiled.exact_source_rows[0].denominator[0].is_zero());
            }
        }
    }
    Ok(())
}

#[test]
fn zero_factor_is_rejected_before_nonzero_constants_are_skipped() -> Result<()> {
    let zero = source::Gaussian::new(Rational::zero(), Rational::zero());
    assert!(matches!(
        square_free_pole_polynomial(vec![zero], variable()),
        Err(Error::Numerical(_))
    ));
    let content = source::Gaussian::new(Rational::from(7), Rational::from(-3));
    assert!(square_free_pole_polynomial(vec![content], variable())?.is_one());
    Ok(())
}

#[test]
fn exact_specialization_uses_stored_dyadics_independently_of_working_profile() -> Result<()> {
    let storage = Precision { bits: 700 };
    let value = storage.rational(&Rational::from((1, 3)));
    let stored = value.re.to_rational();
    assert_ne!(stored, Rational::from((1, 3)));
    let t = Atom::var(variable());
    let a = parameter("repeated_poles_tests::stored_a");
    // Expanded extraction must compute a*a exactly, before numerical rounding.
    let factor = t.pow(2) - Atom::num(2) * &a * &t + a.clone().pow(2);
    let values = [(a, value.clone())].into_iter().collect();
    let expected = exact_polynomial(vec![
        source::Gaussian::new(-stored, Rational::zero()),
        source::Gaussian::new(Rational::one(), Rational::zero()),
    ]);
    let context = RunContext::default();
    let prepared = PreparedRows::new(variable(), &[vec![Atom::one() / &factor]], &context)?;
    for bits in [80, 201, 482] {
        let p = Precision { bits };
        assert_eq!(square_free(&factor, p, &values)?, expected);
        let compiled = prepared.compile(p, &values, &context)?;
        assert_eq!(compiled.poles, vec![p.round(&value)]);
        assert_eq!(compiled.source.values, values);
    }
    Ok(())
}

#[test]
fn tiny_nonzero_dyadic_perturbation_is_coprime_even_if_numeric_arrays_agree() -> Result<()> {
    let storage = Precision { bits: 900 };
    let t = Atom::var(variable());
    let a = parameter("repeated_poles_tests::perturbed_a");
    let d = parameter("repeated_poles_tests::perturbed_d");
    let repeated = t.pow(2) - Atom::num(2) * &a * &t + a.clone().pow(2);
    let perturbed = &repeated + &d;
    let delta = storage.powi(&storage.i(2), -700);
    assert_ne!(delta, storage.zero());
    let values = [(a, storage.rational(&Rational::from((1, 3)))), (d, delta)]
        .into_iter()
        .collect();
    for bits in [80, 482] {
        let p = Precision { bits };
        let exact = specialize(&perturbed, p, &values)?;
        assert!(exact.gcd_euclidean(&exact.derivative()).is_one());
        assert_eq!(square_free(&perturbed, p, &values)?, exact.make_monic());
        assert_eq!(square_free(&perturbed, p, &values)?.degree(), 2);
        assert_eq!(square_free(&repeated, p, &values)?.degree(), 1);
        assert_eq!(
            polynomial_coefficients(&perturbed, variable(), p, &values)?,
            polynomial_coefficients(&repeated, variable(), p, &values)?
        );
    }
    Ok(())
}

#[test]
fn original_rows_factors_content_and_enclosures_are_preserved() -> Result<()> {
    let p = Precision { bits: 482 };
    let context = RunContext::default();
    let values = Default::default();
    let exact_values = source::ExactSpecialization::new(p, &values)?;
    for (denominator, _, signature) in quadratic_cases() {
        let prepared = PreparedRows::new(variable(), &[vec![Atom::one() / denominator]], &context)?;
        let snapshot = format!("{prepared:?}");
        let row = &prepared.rows[0];
        let rational = &row.matrix[0];
        let expected_factors = row
            .denominator_factors
            .iter()
            .map(|(factor, k)| Ok((exact_values.polynomial(factor, variable())?, *k)))
            .collect::<Result<Vec<_>>>()?;
        let compiled = prepared.compile(p, &values, &context)?;
        assert_eq!(format!("{prepared:?}"), snapshot);
        assert_eq!(
            compiled.matrix[0][0].numerator,
            polynomial_coefficients(&rational.numerator, variable(), p, &values)?
        );
        assert_eq!(
            compiled.matrix[0][0].denominator,
            polynomial_coefficients(&rational.denominator, variable(), p, &values)?
        );
        assert_eq!(
            compiled.polynomial_rows[0].denominator,
            polynomial_coefficients(&row.denominator, variable(), p, &values)?
        );
        assert_eq!(
            compiled.exact_source_rows[0].denominator,
            exact_values.polynomial(&row.denominator, variable())?
        );
        assert_eq!(
            compiled.exact_source_rows[0].denominator_factors,
            expected_factors
        );
        let expected_poles = rational
            .new_pole_factors
            .iter()
            .filter(|f| exact_values.polynomial(f, variable()).unwrap().len() > 1)
            .cloned()
            .collect::<Vec<_>>();
        assert_eq!(compiled.pole_polynomials, expected_poles);
        let signatures = expected_factors
            .iter()
            .filter(|(f, _)| f.len() > 1)
            .map(|(f, k)| (f.len() - 1, *k))
            .collect::<Vec<_>>();
        assert_eq!(signatures, vec![signature]);
        let mut product = exact_polynomial(vec![source::Gaussian::new(
            Rational::one(),
            Rational::zero(),
        )]);
        for (factor, multiplicity) in &expected_factors {
            product = &product * &exact_polynomial(factor.clone()).pow(*multiplicity);
        }
        assert_eq!(
            product,
            exact_polynomial(compiled.exact_source_rows[0].denominator.clone())
        );
        let residual = source::ExactSourceResidual::new(
            p,
            &compiled.exact_source_rows,
            &p.zero(),
            &[vec![p.i(1)]],
            1,
        )?;
        let enclosure = &residual.denominators[0];
        let ball = |coefficients: &[source::Gaussian]| {
            source::BallPolynomial::from_coefficients(
                &FloatField::from_rep(source::dyadic_ball(p, &p.zero())),
                coefficients
                    .iter()
                    .map(|c| source::exact_ball(c, p))
                    .collect(),
                Arc::new(PolyVariable::Temporary(0)),
            )
        };
        assert_eq!(
            enclosure.expanded,
            ball(&compiled.exact_source_rows[0].denominator)
        );
        assert_eq!(
            enclosure.factors,
            expected_factors
                .iter()
                .map(|(f, k)| (ball(f), *k))
                .collect::<Vec<_>>()
        );
        let mut product = enclosure.expanded.one();
        for (factor, multiplicity) in &enclosure.factors {
            product = &product * &factor.pow(*multiplicity);
        }
        assert_eq!(product, enclosure.expanded);
    }
    Ok(())
}

#[test]
fn declared_zero_source_guard_survives_numerator_cancellation() -> Result<()> {
    let p = Precision { bits: 482 };
    let j = Atom::var(crate::family::imaginary_parameter());
    let guard = j.pow(2) + 1;
    let mut compiled = DifferentialSystem {
        variable: variable(),
        matrix: vec![vec![&guard / &guard]],
    }
    .compile(p, &Default::default())?;
    assert_eq!(compiled.matrix[0][0].denominator, vec![p.i(1)]);
    assert!(matches!(
        compiled.exclude_polynomials(&[guard]),
        Err(Error::Numerical(_))
    ));
    // A nonzero repeated guard retains its original polynomial, not S.
    let (guard, _, _) = quadratic_cases().into_iter().nth(2).unwrap();
    compiled.exclude_polynomials(&[guard.clone()])?;
    assert_eq!(compiled.poles, vec![p.complex(1, -1)]);
    assert_eq!(compiled.pole_polynomials.len(), 1);
    assert_eq!(
        formal_rational(&compiled.pole_polynomials[0])?,
        formal_rational(&guard)?
    );
    Ok(())
}

#[test]
fn specialized_identical_factors_share_a_pole_but_nearby_distinct_factors_fail_closed() -> Result<()>
{
    let p = Precision { bits: 482 };
    let t = Atom::var(variable());
    let a = parameter("repeated_poles_tests::nearby_a");
    let b = parameter("repeated_poles_tests::nearby_b");
    let prepared = PreparedRows::new(
        variable(),
        &[vec![Atom::one() / (&t - &a), Atom::one() / (&t - &b)]],
        &RunContext::default(),
    )?;
    for (delta_bits, expected_poles) in [(None, Some(1)), (Some(60), Some(2)), (Some(360), None)] {
        let delta = delta_bits
            .map(|b| p.powi(&p.i(2), -b))
            .unwrap_or_else(|| p.zero());
        let values = [(a.clone(), p.i(1)), (b.clone(), p.add(&p.i(1), &delta))]
            .into_iter()
            .collect();
        let result = prepared.compile(p, &values, &RunContext::default());
        if let Some(count) = expected_poles {
            assert_eq!(result?.poles.len(), count);
        } else {
            // Precision::tolerance uses decimal digits: 2^-360 is nonzero
            // at 482 bits, but smaller than the collector's 10^-96 threshold.
            assert_ne!(delta, p.zero());
            assert!(p.norm(&delta) < p.tolerance(p.bits / 5));
            assert!(matches!(result, Err(Error::Numerical(message))
                if message.contains("ambiguous nearby poles")));
        }
    }
    Ok(())
}

#[test]
fn partially_shared_roots_of_unequal_factors_preserve_the_union() -> Result<()> {
    let p = Precision { bits: 482 };
    let t = Atom::var(variable());
    let a = parameter("repeated_poles_tests::shared_a");
    // The quadratic is one formal factor until a=0, where it shares t=1
    // with the first factor. Only the remaining linear factor needs solving.
    let prepared = PreparedRows::new(
        variable(),
        &[vec![
            Atom::one() / (&t - 1),
            Atom::one() / (t.pow(2) - Atom::num(3) * &t + 2 + &a),
        ]],
        &RunContext::default(),
    )?;
    let values = [(a, p.zero())].into_iter().collect();
    let compiled = prepared.compile(p, &values, &RunContext::default())?;
    assert_eq!(compiled.poles, vec![p.i(1), p.i(2)]);
    assert_eq!(compiled.pole_polynomials.len(), 2);
    // The row LCM retains the shared root's original multiplicity.
    assert_eq!(compiled.exact_source_rows[0].denominator.len(), 4);
    Ok(())
}

#[test]
fn actual_quadratics_and_quartic_have_the_same_four_pole_union_in_both_scan_orders() -> Result<()> {
    let p = Precision { bits: 482 };
    let cases = quadratic_cases();
    let qb = cases[0].0.clone();
    let qc = cases[2].0.clone();
    let q4 = quartic();
    let expected_roots = [
        p.scale(&p.complex(1, -1), 1, 2),
        p.scale(&p.complex(5, -5), 1, 2),
        p.complex(1, -1),
        p.scale(&p.complex(5, -5), 1, 8),
    ];
    for (denominators, expected_degrees) in [
        ([qb.clone(), qc.clone(), q4.clone()], [1, 1, 2]),
        ([q4.clone(), qb, qc], [4, 0, 0]),
    ] {
        let matrix = denominators
            .iter()
            .enumerate()
            .map(|(i, denominator)| {
                let mut row = vec![Atom::new(); denominators.len()];
                row[i] = Atom::one() / denominator;
                row
            })
            .collect::<Vec<_>>();
        // Ordinary public differential systems only; no algebraic/gauge flow.
        let compiled = DifferentialSystem {
            variable: variable(),
            matrix: matrix.clone(),
        }
        .compile(p, &Default::default())?;
        assert_eq!(compiled.poles.len(), 4);
        for expected in &expected_roots {
            assert_eq!(
                compiled
                    .poles
                    .iter()
                    .filter(|root| p.close(root, expected, 80))
                    .count(),
                1
            );
        }
        assert_eq!(compiled.source.rows, matrix);
        // In the reverse order, B and C have no new roots but their original
        // factor Atoms must still be retained for source-domain recompilation.
        assert_eq!(compiled.pole_polynomials.len(), denominators.len());
        for (actual, expected) in compiled.pole_polynomials.iter().zip(&denominators) {
            assert_eq!(formal_rational(actual)?, formal_rational(expected)?);
        }
        for (row, denominator) in compiled.exact_source_rows.iter().zip(&denominators) {
            assert_eq!(
                row.denominator,
                specialize(denominator, p, &Default::default())?
                    .coefficients()
                    .to_vec()
            );
        }
        let mut disjoint = Vec::new();
        for (denominator, expected_degree) in denominators.iter().zip(expected_degrees) {
            let original = square_free(denominator, p, &Default::default())?;
            let remaining = unscanned_pole_polynomial(&original, &disjoint)?;
            assert_eq!(remaining.degree(), expected_degree);
            assert!(
                disjoint
                    .iter()
                    .all(|(previous, _)| remaining.gcd_euclidean(previous).is_one())
            );
            if !remaining.is_constant() {
                disjoint.push((remaining, Vec::new()));
            }
        }
        let mut product = specialize(&Atom::one(), p, &Default::default())?;
        for (piece, _) in disjoint {
            product = &product * &piece;
        }
        assert_eq!(product, square_free(&q4, p, &Default::default())?);
    }
    Ok(())
}

#[test]
fn exact_coefficient_extraction_keeps_degrees_above_u16() -> Result<()> {
    let p = Precision { bits: 80 };
    let expression = Atom::var(variable()).pow(100_000) + 1;
    let exact = specialize(&expression, p, &Default::default())?;
    assert_eq!(exact.degree(), 100_000);
    assert!(exact.coefficients()[0].is_one());
    assert!(exact.coefficients()[100_000].is_one());
    assert!(
        exact.coefficients()[1..100_000]
            .iter()
            .all(SingleFloat::is_zero)
    );
    Ok(())
}
