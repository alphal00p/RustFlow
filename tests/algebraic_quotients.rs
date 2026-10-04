use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{
    AlgebraicKinematicSystem, AlgebraicSystem, RootSeed, SquareRoot,
};
use symbolica_amflow::{
    BoundaryData, DifferentialSystem, Error, FlowOptions, Precision, Result, RunContext,
};
use symbolica_amflow::{diffexp, kinematics};

#[test]
fn root_sum_ode_matches_analytic_solutions_on_both_sheets() -> Result<()> {
    let x = symbol!("root_sum_test::x");
    let r = symbol!("root_sum_test::r");
    let roots = vec![SquareRoot {
        symbol: r,
        radicand: Atom::var(x),
    }];
    let p = Precision::decimal(60)?;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 40,
        series_order: 64,
        ..Default::default()
    };
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::one() / (Atom::var(r) + 1)]],
        },
        roots.clone(),
    );
    let compiled = system.compile(p)?;
    let boundary = BoundaryData {
        point: p.i(4),
        values: vec![p.i(1)],
    };
    for (seed, expected) in [
        (RootSeed::Principal, p.scale(&p.exp(&p.i(2)), 9, 16)),
        (RootSeed::Opposite, p.scale(&p.exp(&p.i(-2)), 1, 4)),
    ] {
        let answer = compiled.transport_ordinary(
            &boundary,
            &[p.i(9)],
            &BTreeMap::from([(r, seed)]),
            &options,
            &RunContext::default(),
            false,
        )?;
        assert!(p.norm(&p.sub(&answer.solution.coefficients[0][0], &expected)) < p.tolerance(40));
    }
    // The norm pole is retained even after multiplication cancels it from
    // the numerical kernel. A zero-length request must still enforce it.
    let canceled = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![(Atom::var(x) - 1) / (Atom::var(r) + 1)]],
        },
        roots,
    )
    .compile(p)?;
    assert!(
        canceled
            .singularities()
            .iter()
            .any(|v| p.close(v, &p.i(1), 40))
    );
    assert!(
        canceled
            .transport_ordinary(
                &BoundaryData {
                    point: p.i(1),
                    values: vec![p.i(1)]
                },
                &[p.i(1)],
                &BTreeMap::from([(r, RootSeed::Principal)]),
                &options,
                &RunContext::default(),
                false
            )
            .is_err()
    );
    Ok(())
}

#[test]
fn two_root_dlog_matches_exact_letter_ratio_on_four_sheets() -> Result<()> {
    let x = symbol!("root_sum_dlog::x");
    let s = symbol!("root_sum_dlog::s");
    let r = symbol!("root_sum_dlog::r");
    let q = symbol!("root_sum_dlog::q");
    let eps = symbol!("root_sum_dlog::eps");
    let roots = vec![
        SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        },
        SquareRoot {
            symbol: q,
            radicand: Atom::var(s) + 1,
        },
    ];
    let kinematic = AlgebraicKinematicSystem::canonical_dlog(
        eps,
        &[s],
        &[Atom::var(r) + Atom::var(q) + 1],
        &[vec![vec![Atom::one()]]],
        roots,
    )?;
    let path = kinematics::KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(x))]),
    };
    let system = kinematic.pullback(&path, 3)?;
    for working in [45, 65] {
        let p = Precision::decimal(working)?;
        let compiled = system.compile(p)?;
        let options = FlowOptions {
            digits: 20,
            guard_digits: working - 20,
            series_order: 64,
            ..Default::default()
        };
        let boundary = diffexp::EpsilonBoundary {
            point: p.i(4),
            leading: 0,
            coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
        };
        for sr in [-1, 1] {
            for sq in [-1, 1] {
                let seed = |sign| {
                    if sign == 1 {
                        RootSeed::Principal
                    } else {
                        RootSeed::Opposite
                    }
                };
                let answer = compiled.transport(
                    &boundary,
                    &[p.i(9)],
                    &BTreeMap::from([(r, seed(sr)), (q, seed(sq))]),
                    &options,
                    &RunContext::default(),
                    false,
                )?;
                let a = p.add(
                    &p.i(1 + 2 * sr),
                    &p.scale(&p.pow(&p.i(5), &p.scale(&p.i(1), 1, 2)), sq, 1),
                );
                let b = p.add(
                    &p.i(1 + 3 * sr),
                    &p.scale(&p.pow(&p.i(10), &p.scale(&p.i(1), 1, 2)), sq, 1),
                );
                let log_ratio = p.log(&p.div(&b, &a));
                let mut expected = p.i(1);
                for order in 0..4 {
                    assert!(
                        p.norm(&p.sub(&answer.solution.coefficients[order][0], &expected))
                            < p.tolerance(28)
                    );
                    expected = p.scale(&p.mul(&expected, &log_ratio), 1, order as i64 + 1);
                }
            }
        }
    }
    Ok(())
}

#[test]
fn formal_nonunits_are_rejected_before_stationary_cancellation() -> Result<()> {
    let x = symbol!("root_nonunit::x");
    let t = symbol!("root_nonunit::t");
    let s = symbol!("root_nonunit::s");
    let r = symbol!("root_nonunit::r");
    let eps = symbol!("root_nonunit::eps");
    let system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::one() / (Atom::var(r) - Atom::var(t))]]),
            (t, vec![vec![Atom::new()]]),
        ]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(t).pow(2),
        }],
    };
    let path = kinematics::KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::new()), (t, Atom::var(x))]),
    };
    assert!(matches!(
        system.pullback(&path, 2),
        Err(Error::Unsupported(_))
    ));
    Ok(())
}

#[test]
fn source_domains_survive_stationary_rows_and_radicand_cancellation() -> Result<()> {
    let x = symbol!("root_domain::x");
    let s = symbol!("root_domain::s");
    let t = symbol!("root_domain::t");
    let r = symbol!("root_domain::r");
    let eps = symbol!("root_domain::eps");
    let mut system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::one() / (Atom::var(r) + 1)]]),
            (t, vec![vec![Atom::new()]]),
        ]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(t),
        }],
    };
    let path = kinematics::KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::new()), (t, Atom::var(x))]),
    };
    let p = Precision::decimal(45)?;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 25,
        ..Default::default()
    };
    let reject_at_one = |pulled: AlgebraicSystem| -> Result<()> {
        let compiled = pulled.compile(p)?;
        assert!(
            compiled
                .singularities()
                .iter()
                .any(|v| p.close(v, &p.i(1), 30))
        );
        let boundary = diffexp::EpsilonBoundary {
            point: p.i(1),
            leading: 0,
            coefficients: vec![vec![p.i(1)], vec![p.zero()]],
        };
        assert!(
            compiled
                .transport(
                    &boundary,
                    &[p.i(1)],
                    &BTreeMap::from([(r, RootSeed::Principal)]),
                    &options,
                    &RunContext::default(),
                    false
                )
                .is_err()
        );
        Ok(())
    };
    reject_at_one(system.pullback(&path, 1)?)?;
    // An ordinary rational source pole also survives a zero path Jacobian.
    system.derivatives.get_mut(&s).unwrap()[0][0] = Atom::one() / (Atom::var(t) - 1);
    reject_at_one(system.pullback(&path, 1)?)?;
    // Equal numerator/denominator along the path cannot erase source holes.
    system.derivatives.get_mut(&s).unwrap()[0][0] = Atom::new();
    system.roots[0].radicand = (Atom::var(s) - 1) / (Atom::var(t) - 1);
    let diagonal = kinematics::KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(x)), (t, Atom::var(x))]),
    };
    reject_at_one(system.pullback(&diagonal, 1)?)?;
    Ok(())
}

#[test]
fn mixed_epsilon_norm_keeps_generic_valuation_before_stationary_specialization() -> Result<()> {
    let x = symbol!("mixed_epsilon_domain::x");
    let s = symbol!("mixed_epsilon_domain::s");
    let t = symbol!("mixed_epsilon_domain::t");
    let r = symbol!("mixed_epsilon_domain::r");
    let eps = symbol!("mixed_epsilon_domain::eps");
    let system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([
            (
                s,
                vec![vec![Atom::one() / (Atom::var(r) + 1 + Atom::var(eps))]],
            ),
            (t, vec![vec![Atom::new()]]),
        ]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(t) + 1,
        }],
    };
    let path = |t_value| kinematics::KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::new()), (t, t_value)]),
    };
    // Norm=(1+eps)^2-(t+1)=2eps+eps^2-t. Its generic leading
    // coefficient is -t. The t=0 path must not rebase that valuation to eps^1
    // and discard the source domain merely because the s Jacobian is zero.
    assert!(matches!(
        system.pullback(&path(Atom::new()), 2),
        Err(Error::InvalidInput(_))
    ));
    let moving = system.pullback(&path(Atom::var(x)), 2)?;
    assert!(moving.nonzero_conditions.iter().any(|guard| {
        guard
            .replace(Atom::var(x))
            .with(Atom::new())
            .together()
            .cancel()
            .is_zero()
    }));
    // Both sheets are regular at t=3, and stationary transport remains valid.
    let regular = system.pullback(&path(Atom::num(3)), 2)?;
    assert!(
        regular
            .system
            .matrices
            .iter()
            .flatten()
            .flatten()
            .all(|a| a.is_zero())
    );
    let p = Precision::decimal(45)?;
    let compiled = regular.compile(p)?;
    let options = FlowOptions {
        digits: 20,
        guard_digits: 25,
        ..Default::default()
    };
    let boundary = diffexp::EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]],
    };
    for seed in [RootSeed::Principal, RootSeed::Opposite] {
        let answer = compiled.transport(
            &boundary,
            &[p.i(1)],
            &BTreeMap::from([(r, seed)]),
            &options,
            &RunContext::default(),
            false,
        )?;
        assert_eq!(answer.solution.coefficients, boundary.coefficients);
    }
    Ok(())
}

#[test]
fn public_source_guards_keep_exact_complex_coefficients() -> Result<()> {
    let s = symbol!("complex_norm_domain::s");
    let r = symbol!("complex_norm_domain::r");
    let eps = symbol!("complex_norm_domain::eps");
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let radicand = Atom::one() + imaginary * Atom::var(s) * (Atom::var(s) - 1);
    let system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]])]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: radicand.clone(),
        }],
    };
    let conditions = system.nonzero_conditions()?;
    assert!(
        conditions
            .iter()
            .any(|a| (a - &radicand).together().cancel().is_zero())
    );
    assert!(conditions.iter().all(|a| {
        a.derivative(symbol!("symbolica_amflow::imaginary_unit"))
            .is_zero()
    }));
    Ok(())
}
