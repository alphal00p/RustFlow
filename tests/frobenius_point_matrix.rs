//! Streamed point-matrix regressions with independent exp/log oracles.
use std::sync::Arc;
use symbolica::prelude::*;
use symbolica_amflow::{numeric::solve, *};

fn weighted(
    p: Precision,
    matrix: &[Vec<ComplexFloat>],
    constants: &[ComplexFloat],
) -> Vec<ComplexFloat> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .zip(constants)
                .fold(p.zero(), |sum, (a, b)| p.add(&sum, &p.mul(a, b)))
        })
        .collect()
}

fn branch_log(p: Precision, point: &ComplexFloat, winding: i32) -> ComplexFloat {
    p.add(
        &p.log(point),
        &p.scale(&p.log(&p.i(-1)), 2 * i64::from(winding), 1),
    )
}

fn logarithmic_oracle(
    p: Precision,
    point: &ComplexFloat,
    winding: i32,
    a: &ComplexFloat,
    b: &ComplexFloat,
    c: Option<&ComplexFloat>,
) -> Vec<ComplexFloat> {
    let log = branch_log(p, point, winding);
    let power = p.exp(&p.scale(&log, 1, 5));
    let mut values = vec![p.mul(&power, &p.add(a, &p.mul(b, &log))), p.mul(&power, b)];
    if let Some(c) = c {
        let inner = p.add(
            c,
            &p.add(
                &p.mul(a, &log),
                &p.scale(&p.mul(b, &p.mul(&log, &log)), 1, 2),
            ),
        );
        values.push(p.mul(&p.mul(point, &power), &inner));
    }
    values
}

fn assert_relative(p: Precision, actual: &ComplexFloat, expected: &ComplexFloat, digits: u32) {
    assert!(p.finite(actual) && p.finite(expected));
    if *expected == p.zero() {
        assert!(p.close(actual, expected, digits));
    } else {
        assert!(p.norm(&p.sub(actual, expected)) / p.norm(expected) < p.tolerance(digits));
    }
}

#[test]
fn streamed_point_matrix_preserves_columns_gauge_winding_and_multiple_rhs() -> Result<()> {
    let x = symbol!("point_matrix_columns::x");
    let z = Atom::var(x);
    let lambda = Atom::num(Rational::from((1, 5)));
    let nonresonant = vec![
        vec![&lambda / &z, Atom::one() / &z],
        vec![Atom::new(), &lambda / &z],
    ];
    let resonant = vec![
        vec![&lambda / &z, Atom::one() / &z, Atom::new()],
        vec![Atom::new(), &lambda / &z, Atom::new()],
        vec![Atom::one(), Atom::new(), (&lambda + 1) / &z],
    ];
    let sheared = vec![
        vec![z.clone().pow(-2), z.clone().pow(-2)],
        vec![-z.clone().pow(-2), -z.clone().pow(-2)],
    ];
    for (fixture, matrix) in [nonresonant, resonant, sheared].into_iter().enumerate() {
        let prepared = DifferentialSystem {
            variable: x,
            matrix,
        }
        .prepare_frobenius(&RunContext::default())?;
        for digits in [80_u32, 120] {
            let p = Precision::decimal(digits)?;
            let parameters = Default::default();
            let full = prepared.evaluate(p, &parameters, 32, &RunContext::default())?;
            if fixture < 2 {
                assert!(
                    full.columns
                        .iter()
                        .any(|column| column.coefficients.iter().any(|logs| logs.len() > 1))
                );
            } else {
                assert!(
                    prepared
                        .basis_transformation()
                        .iter()
                        .enumerate()
                        .any(|(i, row)| {
                            row.iter().enumerate().any(|(j, a)| {
                                !(a.clone() - Atom::num(i64::from(i == j)))
                                    .together()
                                    .cancel()
                                    .is_zero()
                            })
                        })
                );
            }
            for point in [
                p.rational(&Rational::from((1, 8))),
                p.parse("0.05", "0.07")?,
            ] {
                for winding in [-1_i32, 0, 2] {
                    // Generate the streamed matrix once, then reuse it for four RHS.
                    let small = prepared.evaluate_matrix_with_winding(
                        p,
                        &parameters,
                        32,
                        &point,
                        winding,
                        &RunContext::default(),
                    )?;
                    let legacy = full.evaluate_with_winding(&point, &parameters, winding)?;
                    assert_eq!(small.len(), legacy.len());
                    for (row, expected) in small.iter().zip(&legacy) {
                        assert_eq!(row.len(), expected.len());
                        for (actual, expected) in row.iter().zip(expected) {
                            assert!(p.close(actual, expected, digits - 20));
                        }
                    }
                    for rhs in 0..4_i64 {
                        let constants = (0..small.len())
                            .map(|i| p.complex(rhs + i as i64 + 2, rhs - i as i64 - 1))
                            .collect::<Vec<_>>();
                        let target = weighted(p, &legacy, &constants);
                        let recovered = solve(p, small.clone(), target.clone())?;
                        // FrobeniusBasis::match_values is principal-branch only.
                        if winding == 0 {
                            let old = full.match_values(&point, &target, &parameters)?;
                            for (a, b) in recovered.iter().zip(old) {
                                assert!(p.close(a, &b, digits - 25));
                            }
                        }
                        let direct = prepared.match_values_with_winding(
                            p,
                            &parameters,
                            32,
                            &point,
                            &target,
                            winding,
                            &RunContext::default(),
                        )?;
                        for ((actual, expected), direct) in
                            recovered.iter().zip(&constants).zip(&direct)
                        {
                            assert!(p.close(actual, expected, digits - 25));
                            assert!(p.close(actual, direct, digits - 25));
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn streamed_point_matrix_matches_independent_exp_log_and_tiny_rhs_oracles() -> Result<()> {
    let x = symbol!("point_matrix_oracle::x");
    let z = Atom::var(x);
    let lambda = Atom::num(Rational::from((1, 5)));
    for resonant in [false, true] {
        let mut matrix = vec![
            vec![&lambda / &z, Atom::one() / &z],
            vec![Atom::new(), &lambda / &z],
        ];
        if resonant {
            matrix[0].push(Atom::new());
            matrix[1].push(Atom::new());
            matrix.push(vec![Atom::one(), Atom::new(), (&lambda + 1) / &z]);
        }
        let prepared = DifferentialSystem {
            variable: x,
            matrix,
        }
        .prepare_frobenius(&RunContext::default())?;
        let p = Precision::decimal(160)?;
        let parameters = Default::default();
        let anchor = p.rational(&Rational::from((1, 8)));
        let endpoint = p.parse("1e-10", "0")?;
        for winding in [-1_i32, 0, 2] {
            let at_anchor = prepared.evaluate_matrix_with_winding(
                p,
                &parameters,
                24,
                &anchor,
                winding,
                &RunContext::default(),
            )?;
            let at_endpoint = prepared.evaluate_matrix_with_winding(
                p,
                &parameters,
                24,
                &endpoint,
                winding,
                &RunContext::default(),
            )?;
            for tiny in [false, true] {
                let a = if tiny { p.zero() } else { p.complex(2, -1) };
                let b = if tiny {
                    p.parse("1e-120", "0")?
                } else {
                    p.complex(-3, 2)
                };
                let c = if tiny { p.zero() } else { p.complex(1, 1) };
                let target =
                    logarithmic_oracle(p, &anchor, winding, &a, &b, resonant.then_some(&c));
                let constants = solve(p, at_anchor.clone(), target)?;
                let actual = weighted(p, &at_endpoint, &constants);
                let expected =
                    logarithmic_oracle(p, &endpoint, winding, &a, &b, resonant.then_some(&c));
                for (a, b) in actual.iter().zip(&expected) {
                    assert_relative(p, a, b, 80);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn streamed_point_matrix_keeps_domain_order_precision_and_cancellation_guards() -> Result<()> {
    let x = symbol!("point_matrix_guards::x");
    let e = parse!("point_matrix_guards::epsilon");
    let prepared = DifferentialSystem {
        variable: x,
        matrix: vec![vec![parse!("(point_matrix_guards::epsilon^2-1)/((point_matrix_guards::epsilon-1)*point_matrix_guards::x)")]],
    }.prepare_frobenius(&RunContext::default())?;
    let p = Precision::decimal(80)?;
    let point = p.rational(&Rational::from((1, 2)));
    let excluded = ahash::HashMap::from_iter([(e.clone(), p.i(1))]);
    assert!(matches!(
        prepared.evaluate_matrix_with_winding(p, &excluded, 8, &point, 0, &RunContext::default()),
        Err(Error::Numerical(_))
    ));
    let parameters = ahash::HashMap::from_iter([(e, p.i(2))]);
    assert!(matches!(
        prepared.evaluate_matrix_with_winding(
            p,
            &parameters,
            8,
            &p.zero(),
            0,
            &RunContext::default()
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        prepared.evaluate_matrix_with_winding(
            Precision { bits: u32::MAX },
            &parameters,
            8,
            &point,
            0,
            &RunContext::default()
        ),
        Err(Error::Limit(_))
    ));
    assert!(matches!(
        prepared.evaluate_matrix_with_winding(
            p,
            &parameters,
            200_001,
            &point,
            0,
            &RunContext::default()
        ),
        Err(Error::Limit(_))
    ));
    assert!(matches!(
        prepared.match_values(p, &parameters, 8, &point, &[], &RunContext::default()),
        Err(Error::InvalidInput(_))
    ));
    let token = CancellationToken::default();
    let cancel = token.clone();
    let context = RunContext {
        cancellation: token,
        progress: Some(Arc::new(move |_| cancel.cancel())),
    };
    assert!(matches!(
        prepared.evaluate_matrix_with_winding(p, &parameters, 8, &point, 0, &context),
        Err(Error::Cancelled)
    ));
    let values = prepared.evaluate_matrix_with_winding(
        p,
        &parameters,
        8,
        &point,
        0,
        &RunContext::default(),
    )?;
    assert_eq!(values.len(), 1);
    assert_eq!(values[0].len(), 1);
    assert!(p.finite(&values[0][0]));

    let eps = parse!("point_matrix_guards::near_epsilon");
    let z = Atom::var(x);
    let near = DifferentialSystem {
        variable: x,
        matrix: vec![
            vec![&eps / &z, Atom::new(), Atom::new()],
            vec![Atom::new(), &eps / &z, Atom::new()],
            vec![Atom::new(), Atom::new(), Atom::num(2) / &z],
        ],
    }
    .prepare_frobenius(&RunContext::default())?;
    let parameters = ahash::HashMap::from_iter([(eps, p.parse("1e-100", "0")?)]);
    assert!(matches!(
        near.evaluate_matrix_with_winding(p, &parameters, 4, &point, 0, &RunContext::default()),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}
