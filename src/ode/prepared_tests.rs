//! Frozen legacy compiler oracle plus independent source-specialization checks.
use super::*;

fn same(actual: &CompiledSystem, legacy: &CompiledSystem) {
    assert_eq!(actual.p.bits, legacy.p.bits);
    assert_eq!(actual.source.variable, legacy.source.variable);
    assert_eq!(actual.source.rows, legacy.source.rows);
    assert_eq!(actual.source.values, legacy.source.values);
    assert_eq!(actual.poles, legacy.poles);
    assert_eq!(actual.pole_polynomials, legacy.pole_polynomials);
    for (a, b) in actual
        .matrix
        .iter()
        .flatten()
        .zip(legacy.matrix.iter().flatten())
    {
        assert_eq!(a.numerator, b.numerator);
        assert_eq!(a.denominator, b.denominator);
        for value in a.numerator.iter().chain(&a.denominator) {
            assert_eq!(value.re.as_raw().prec(), actual.p.bits);
            assert_eq!(value.im.as_raw().prec(), actual.p.bits);
        }
    }
    for (a, b) in actual.polynomial_rows.iter().zip(&legacy.polynomial_rows) {
        assert_eq!(a.denominator, b.denominator);
        assert_eq!(a.entries, b.entries);
    }
    for (a, b) in actual
        .exact_source_rows
        .iter()
        .zip(&legacy.exact_source_rows)
    {
        assert_eq!(a.denominator, b.denominator);
        assert_eq!(a.denominator_factors, b.denominator_factors);
        assert_eq!(a.entries, b.entries);
    }
}

#[test]
fn prepared_rows_preserve_legacy_complex_coefficients_factors_poles_and_source() {
    let x = symbol!("prepared_rows_tests::x");
    let a = symbol!("prepared_rows_tests::a");
    let x_atom = Atom::var(x);
    let parameter = Atom::var(a);
    let i = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let rows = vec![
        vec![
            (&parameter + &i * &x_atom) / (Atom::num(6) * (&x_atom - 1).pow(2)),
            (&parameter - &x_atom) / ((&x_atom - 1) * (&x_atom + 2)),
            Atom::new(),
        ],
        vec![Atom::num(3) / (Atom::num(7) * (&x_atom + &i)), Atom::new()],
        vec![Atom::new()],
    ];
    let context = RunContext::default();
    let prepared = prepared::PreparedRows::new(x, &rows, &context).unwrap();
    let storage = Precision { bits: 700 };
    for value in [
        storage.rational(&Rational::from((1, 3))),
        storage.complex(2, -3),
        storage.zero(),
    ] {
        let values = [(parameter.clone(), value)].into_iter().collect();
        for bits in [80, 201, 400] {
            let p = Precision { bits };
            let actual = prepared.compile(p, &values, &context).unwrap();
            let legacy = legacy_compile_rows(x, &rows, p, &values).unwrap();
            same(&actual, &legacy);
            assert!(actual.poles.iter().any(|z| p.close(z, &p.i(1), 15)));
            assert!(actual.poles.iter().any(|z| p.close(z, &p.i(-2), 15)));
            assert!(
                actual
                    .poles
                    .iter()
                    .any(|z| p.close(z, &p.complex(0, -1), 15))
            );
        }
    }
}

#[test]
fn prepared_rows_snapshot_and_stored_parameter_precision_are_independent() {
    let x = symbol!("prepared_rows_tests::snapshot_x");
    let a = Atom::var(symbol!("prepared_rows_tests::snapshot_a"));
    let mut rows = vec![vec![a.clone()]];
    let context = RunContext::default();
    let prepared = prepared::PreparedRows::new(x, &rows, &context).unwrap();
    rows[0][0] = Atom::num(99);
    let storage = Precision { bits: 700 };
    let value = storage.rational(&Rational::from((1, 3)));
    let low = Precision { bits: 80 };
    let high = Precision { bits: 400 };
    let values = [(a.clone(), value.clone())].into_iter().collect();
    let first = prepared.compile(low, &values, &context).unwrap();
    let second = prepared.compile(high, &values, &context).unwrap();
    assert_eq!(first.source.rows, vec![vec![a.clone()]]);
    assert_eq!(second.matrix[0][0].numerator[0], high.round(&value));
    assert_ne!(
        second.matrix[0][0].numerator[0],
        high.round(&first.matrix[0][0].numerator[0])
    );
    assert_eq!(
        first.exact_source_rows[0].entries[0].1[0].re,
        value.re.to_rational()
    );
    assert_eq!(
        second.exact_source_rows[0].entries[0].1[0].re,
        value.re.to_rational()
    );
    let changed = [(a, storage.zero())].into_iter().collect();
    let zero = prepared.compile(high, &changed, &context).unwrap();
    assert!(zero.polynomial_rows[0].entries.is_empty());
    assert!(zero.exact_source_rows[0].entries.is_empty());
    assert_eq!(first.matrix[0][0].numerator[0], low.round(&value));
}

#[test]
fn prepared_rows_preserve_parameter_domain_failures_and_cancellation() {
    let x = symbol!("prepared_rows_tests::domain_x");
    let a = Atom::var(symbol!("prepared_rows_tests::domain_a"));
    let rows = vec![vec![Atom::one() / (&a * (Atom::var(x) - 1))]];
    let context = RunContext::default();
    let prepared = prepared::PreparedRows::new(x, &rows, &context).unwrap();
    let p = Precision { bits: 120 };
    let values = [(a, p.zero())].into_iter().collect();
    assert!(matches!(
        prepared.compile(p, &values, &context),
        Err(Error::Numerical(_))
    ));
    assert!(matches!(
        legacy_compile_rows(x, &rows, p, &values),
        Err(Error::Numerical(_))
    ));
    context.cancellation.cancel();
    assert!(matches!(
        prepared.compile(p, &values, &context),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        prepared::PreparedRows::new(x, &rows, &context),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        prepared::PreparedRows::new(x, &[vec![]], &RunContext::default()),
        Err(Error::InvalidInput(_))
    ));
}

// Original production compiler from commit 20a2af0, retained as an independent
// regression oracle while exact preparation is shared across precision profiles.
fn legacy_compile_rows(
    variable: Symbol,
    rows: &[Vec<Atom>],
    p: Precision,
    values: &ahash::HashMap<Atom, C>,
) -> Result<CompiledSystem> {
    validate_ode_variable(variable)?;
    let mut matrix = Vec::new();
    let mut poles = Vec::new();
    let mut pole_polynomials = Vec::new();
    let mut seen_factors = ahash::HashSet::default();
    let mut factorizations = ahash::HashMap::default();
    let mut polynomial_rows = Vec::new();
    let mut exact_source_rows = Vec::new();
    let exact_values = source::ExactSpecialization::new(p, values)?;
    for row in rows {
        let mut out = Vec::new();
        let mut exact = Vec::new();
        let mut common_denominator: Option<ExactPolynomial> = None;
        for a in row {
            let rational: RationalPolynomial<IntegerRing, u16> = crate::family::encode_complex(a)
                .try_to_rational_polynomial(&Q, &Z, None)
                .map_err(|e| Error::InvalidInput(e.to_string()))?;
            let numerator =
                polynomial_coefficients(&rational.numerator.to_expression(), variable, p, values)?;
            let denominator = polynomial_coefficients(
                &rational.denominator.to_expression(),
                variable,
                p,
                values,
            )?;
            if denominator.iter().all(|v| *v == p.zero()) {
                return Err(Error::Numerical(
                    "identically zero specialized denominator".into(),
                ));
            }
            let factors = denominator_factors(&rational.denominator, &mut factorizations)?;
            for (factor, _) in factors.iter() {
                let base = factor.to_expression();
                if !seen_factors.insert(base.clone()) {
                    continue;
                }
                let coefficients = polynomial_coefficients(&base, variable, p, values)?;
                if coefficients.len() <= 1 {
                    continue;
                }
                pole_polynomials.push(base);
                let roots = polynomial_roots(p, &coefficients, variable)?;
                for root in roots {
                    retain_distinct_pole(p, &mut poles, root);
                }
            }
            out.push(NumericRational {
                numerator,
                denominator,
            });
            common_denominator = Some(if let Some(previous) = common_denominator {
                let quotient = previous
                    .try_div(&previous.gcd(&rational.denominator))
                    .ok_or_else(|| {
                        Error::Numerical("exact denominator LCM division failed".into())
                    })?;
                multiply_polynomials(&quotient, &rational.denominator)?
            } else {
                rational.denominator.clone()
            });
            exact.push(rational);
        }
        let common_denominator = common_denominator.unwrap();
        let denominator =
            polynomial_coefficients(&common_denominator.to_expression(), variable, p, values)?;
        let exact_denominator =
            exact_values.polynomial(&common_denominator.to_expression(), variable)?;
        let exact_denominator_factors =
            denominator_factors(&common_denominator, &mut factorizations)?
                .iter()
                .map(|(factor, multiplicity)| {
                    Ok((
                        exact_values.polynomial(&factor.to_expression(), variable)?,
                        *multiplicity,
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
        if exact_denominator_factors
            .iter()
            .any(|(factor, _)| factor.iter().all(SingleFloat::is_zero))
        {
            return Err(Error::Numerical(
                "identically zero specialized denominator factor".into(),
            ));
        }
        let mut exact_entries = Vec::new();
        let mut entries = Vec::new();
        for (j, rational) in exact.iter().enumerate() {
            if rational.numerator.is_zero() {
                continue;
            }
            let multiplier = common_denominator
                .try_div(&rational.denominator)
                .ok_or_else(|| Error::Numerical("exact denominator clearing failed".into()))?;
            let cleared = multiply_polynomials(&rational.numerator, &multiplier)?;
            let coefficients =
                polynomial_coefficients(&cleared.to_expression(), variable, p, values)?;
            let exact_coefficients = exact_values.polynomial(&cleared.to_expression(), variable)?;
            if exact_coefficients.iter().any(|c| !c.is_zero()) {
                exact_entries.push((j, exact_coefficients));
            }
            if coefficients.iter().any(|c| *c != p.zero()) {
                entries.push((j, coefficients));
            }
        }
        polynomial_rows.push(PolynomialRow {
            denominator,
            entries,
        });
        exact_source_rows.push(source::ExactPolynomialRow {
            denominator: exact_denominator,
            denominator_factors: exact_denominator_factors,
            entries: exact_entries,
        });
        matrix.push(out);
    }
    Ok(CompiledSystem {
        source: Arc::new(ExactRows {
            variable,
            rows: rows.to_vec(),
            values: values.clone(),
        }),
        p,
        matrix,
        poles,
        pole_polynomials,
        polynomial_rows,
        exact_source_rows,
    })
}
