#![allow(clippy::needless_range_loop)] // Explicit epsilon/master indices define the exact augmented system.
use std::collections::BTreeMap;
use std::time::Instant;
use symbolica_amflow::AsymptoticConstraint;
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem};
use symbolica_amflow::symbolica::prelude::*;
use symbolica_amflow::{
    ComplexFloat as C, DifferentialSystem, Error, FlowOptions, Precision, Result, RunContext,
};

pub(crate) fn atom(s: &str) -> Result<Atom> {
    Atom::parse(
        s.replace('[', "(").replace(']', ")"),
        "banana_proto",
        Default::default(),
    )
    .map_err(|e| Error::InvalidInput(e.to_string()))
}
pub(crate) fn evaluate_equal(
    digits: u32,
    order: usize,
    endpoint: Rational,
) -> Result<(EpsilonBoundary, EpsilonBoundary)> {
    let fixture: serde_json::Value =
        serde_json::from_str(include_str!("../../fixtures/diffexp/banana-equal.json")).unwrap();
    let p = Precision::decimal(digits)?;
    let physical_matrices = fixture["matrices"]
        .as_array()
        .unwrap()
        .iter()
        .map(|m| {
            m["matrix"]
                .as_array()
                .unwrap()
                .iter()
                .map(|row| {
                    row.as_array()
                        .unwrap()
                        .iter()
                        .map(|x| atom(x.as_str().unwrap()))
                        .collect::<Result<Vec<_>>>()
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let n = 20;
    let mut augmented = vec![vec![Atom::new(); n]; n];
    for k in 0..5 {
        for l in 0..=k {
            for i in 0..4 {
                for j in 0..4 {
                    augmented[4 * k + i][4 * l + j] = physical_matrices[k - l][i][j].clone();
                }
            }
        }
    }
    let t = symbol!("banana_proto::t");
    let x = symbol!("banana_proto::x");
    let physical = DifferentialSystem {
        variable: t,
        matrix: augmented,
    };
    let mapping = BTreeMap::from([(Atom::var(t), -Atom::var(x).pow(-1))]);
    let infinity = DifferentialSystem {
        variable: x,
        matrix: physical
            .matrix
            .iter()
            .map(|r| {
                r.iter()
                    .map(|a| {
                        (symbolica_amflow::family::substitute(a, &mapping) / Atom::var(x).pow(2))
                            .together()
                            .cancel()
                    })
                    .collect()
            })
            .collect(),
    };
    eprintln!("start Frobenius, digits={digits}, order={order}");
    let clock = Instant::now();
    let basis = infinity.frobenius(p, &Default::default(), order)?;
    let frobenius_seconds = clock.elapsed().as_secs_f64();
    eprintln!(
        "Frobenius built in {frobenius_seconds:.3}s, {} columns",
        basis.columns.len()
    );
    let max_logs = basis
        .columns
        .iter()
        .flat_map(|c| c.coefficients.iter())
        .map(Vec::len)
        .max()
        .unwrap();
    let logx = symbol!("banana_proto::logx");
    let mut parameters = ahash::HashMap::default();
    parameters.insert(atom("Pi")?, C::new(p.log(&p.i(-1)).im, p.real(0)));
    for k in [3, 5] {
        parameters.insert(
            atom(&format!("Zeta({k})"))?,
            C::new(
                Float::with_val(p.bits, p.real(k).as_raw().clone().zeta()),
                p.real(0),
            ),
        );
    }
    let leading = fixture["leading_epsilon_coefficients_by_master_wolfram"]
        .as_array()
        .unwrap();
    let leading = leading
        .iter()
        .map(|r| {
            r.as_array()
                .unwrap()
                .iter()
                .map(|x| atom(x.as_str().unwrap()))
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let mut constraints = Vec::new();
    for k in 0..5 {
        for log in 0..max_logs {
            let expansion = leading[2][k]
                .series(logx, 0, max_logs as i64 + 1)
                .map_err(|e| Error::InvalidInput(e.to_string()))?;
            let value = p.eval(
                &expansion
                    .coefficient(Rational::from(log as i64))
                    .unwrap_or_default(),
                &parameters,
            )?;
            constraints.push(AsymptoticConstraint {
                component: 4 * k + 2,
                power: Atom::num(1),
                log_power: log,
                coefficient: value,
            });
            constraints.push(AsymptoticConstraint {
                component: 4 * k + 2,
                power: Atom::num(0),
                log_power: log,
                coefficient: p.zero(),
            });
            constraints.push(AsymptoticConstraint {
                component: 4 * k + 3,
                power: Atom::num(0),
                log_power: log,
                coefficient: if log == 0 {
                    p.eval(&leading[3][k], &parameters)?
                } else {
                    p.zero()
                },
            });
        }
    }
    let matching = Instant::now();
    let constants = basis.match_constraints(&constraints)?;
    let matching_seconds = matching.elapsed().as_secs_f64();
    eprintln!(
        "matched {} constraints in {matching_seconds:.3}s, max_logs={max_logs}",
        constraints.len()
    );
    let z = p.rational(&Rational::from((1, 128)));
    let f = basis.evaluate(&z, &Default::default())?;
    let values: Vec<C> = f
        .iter()
        .map(|row| {
            row.iter()
                .zip(&constants)
                .fold(p.zero(), |s, (a, c)| p.add(&s, &p.mul(a, c)))
        })
        .collect();
    let options = FlowOptions {
        digits: 20,
        guard_digits: digits - 20,
        series_order: order.max(40),
        ..Default::default()
    };
    let infinity_epsilon = EpsilonSystem {
        variable: x,
        matrices: physical_matrices
            .iter()
            .map(|matrix| {
                matrix
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|a| {
                                (symbolica_amflow::family::substitute(a, &mapping)
                                    / Atom::var(x).pow(2))
                                .together()
                                .cancel()
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect(),
    };
    let ordinary = Instant::now();
    let minus1 = infinity_epsilon
        .compile(p, &Default::default())?
        .transport(
            &EpsilonBoundary {
                point: z,
                leading: 0,
                coefficients: values.chunks(4).map(<[C]>::to_vec).collect(),
            },
            &[p.i(1)],
            &options,
            &RunContext::default(),
            false,
        )?;
    let infinity_transport_seconds = ordinary.elapsed().as_secs_f64();
    let compiled = EpsilonSystem {
        variable: t,
        matrices: physical_matrices,
    }
    .compile(p, &Default::default())?;
    let end = p.rational(&endpoint);
    let route = compiled.plan_path(&p.i(-1), &end, 1)?;
    let ordinary = Instant::now();
    let physical_result = compiled.transport(
        &EpsilonBoundary {
            point: p.i(-1),
            leading: 0,
            coefficients: minus1.coefficients.clone(),
        },
        &route,
        &options,
        &RunContext::default(),
        false,
    )?;
    let physical_transport_seconds = ordinary.elapsed().as_secs_f64();
    eprintln!(
        "infinity transport {infinity_transport_seconds:.3}s/{} steps, physical {physical_transport_seconds:.3}s/{} steps",
        minus1.diagnostics.steps, physical_result.diagnostics.steps
    );
    Ok((
        EpsilonBoundary {
            point: p.i(-1),
            leading: 0,
            coefficients: minus1.coefficients,
        },
        EpsilonBoundary {
            point: end,
            leading: 0,
            coefficients: physical_result.coefficients,
        },
    ))
}
