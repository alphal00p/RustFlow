use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, CanonicalAlgebraicSystem, SquareRoot};
use symbolica_amflow::family::substitute;
use symbolica_amflow::kinematics::KinematicPath;
use symbolica_amflow::{Error, Result};

#[test]
fn canonical_path_first_matches_exact_dense_total_derivatives() -> Result<()> {
    let a = symbol!("canonical_equal::a");
    let b = symbol!("canonical_equal::b");
    let r = symbol!("canonical_equal::r");
    let x = symbol!("canonical_equal::x");
    let eps = symbol!("canonical_equal::eps");
    let roots = vec![SquareRoot {
        symbol: r,
        radicand: Atom::var(b).pow(2) - Atom::var(a),
    }];
    let letters = vec![
        (Atom::var(a) + Atom::var(r)) / (Atom::var(a) - Atom::var(r)),
        Atom::var(a) + Atom::var(b),
    ];
    let matrices = vec![
        vec![
            vec![Atom::num(1), Atom::num(2)],
            vec![Atom::new(), Atom::new()],
        ],
        vec![
            vec![Atom::num(-1), Atom::new()],
            vec![Atom::num(3), Atom::num(1)],
        ],
    ];
    let canonical =
        CanonicalAlgebraicSystem::new(eps, &[b, a], &letters, &matrices, roots.clone())?;
    assert_eq!(canonical.dimension(), 2);
    assert_eq!(canonical.variables(), &[b, a]);
    assert_eq!(canonical.epsilon(), eps);
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([
            (a, Atom::num(1) + Atom::var(x)),
            (b, Atom::num(4) + Atom::var(x)),
        ]),
    };
    let direct = canonical.pullback(&path, 3)?;
    let dense = AlgebraicKinematicSystem::canonical_dlog(eps, &[b, a], &letters, &matrices, roots)?
        .pullback(&path, 3)?;
    for (a, b) in direct
        .system
        .matrices
        .iter()
        .flatten()
        .flatten()
        .zip(dense.system.matrices.iter().flatten().flatten())
    {
        assert!((a - b).together().cancel().is_zero());
    }
    assert_eq!(direct.roots[0].symbol, dense.roots[0].symbol);
    assert!(
        (&direct.roots[0].radicand - &dense.roots[0].radicand)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(!canonical.nonzero_conditions().is_empty());
    assert!(!direct.nonzero_conditions.is_empty());
    Ok(())
}

#[test]
fn stationary_jacobians_and_cancelling_letters_keep_source_holes() -> Result<()> {
    let s = symbol!("canonical_stationary::s");
    let t = symbol!("canonical_stationary::t");
    let x = symbol!("canonical_stationary::x");
    let eps = symbol!("canonical_stationary::eps");
    let letter = Atom::var(s) - Atom::var(t);
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s, t],
        &[letter.clone(), letter],
        &[vec![vec![Atom::one()]], vec![vec![Atom::num(-1)]]],
        vec![],
    )?;
    let diagonal = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(x)), (t, Atom::var(x))]),
    };
    assert!(matches!(
        system.pullback(&diagonal, 2),
        Err(Error::InvalidInput(_))
    ));
    let regular = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(x)), (t, Atom::one())]),
    };
    let pulled = system.pullback(&regular, 2)?;
    assert!(
        pulled
            .system
            .matrices
            .iter()
            .flatten()
            .flatten()
            .all(Atom::is_zero)
    );
    let point = BTreeMap::from([(Atom::var(x), Atom::one())]);
    assert!(
        pulled
            .nonzero_conditions
            .iter()
            .any(|a| substitute(a, &point).together().cancel().is_zero())
    );
    let single = CanonicalAlgebraicSystem::new(
        eps,
        &[s, t],
        &[Atom::var(s)],
        &[vec![vec![Atom::one()]]],
        vec![],
    )?;
    let stationary = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::new()), (t, Atom::var(x))]),
    };
    assert!(matches!(
        single.pullback(&stationary, 2),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn radicand_cancellation_does_not_erase_its_original_domain() -> Result<()> {
    let s = symbol!("canonical_radicand::s");
    let t = symbol!("canonical_radicand::t");
    let r = symbol!("canonical_radicand::r");
    let x = symbol!("canonical_radicand::x");
    let eps = symbol!("canonical_radicand::eps");
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s, t],
        &[Atom::one()],
        &[vec![vec![Atom::new()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: (Atom::var(s) - 1) / (Atom::var(t) - 1),
        }],
    )?;
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(x)), (t, Atom::var(x))]),
    };
    let pulled = system.pullback(&path, 1)?;
    assert_eq!(pulled.roots[0].radicand, Atom::one());
    let point = BTreeMap::from([(Atom::var(x), Atom::one())]);
    assert!(
        pulled
            .nonzero_conditions
            .iter()
            .any(|a| substitute(a, &point).together().cancel().is_zero())
    );
    let singular = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::one()), (t, Atom::var(x))]),
    };
    assert!(matches!(
        system.pullback(&singular, 1),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn root_norms_reject_exceptional_paths_even_when_one_sheet_is_regular() -> Result<()> {
    let s = symbol!("canonical_norm::s");
    let t = symbol!("canonical_norm::t");
    let r = symbol!("canonical_norm::r");
    let x = symbol!("canonical_norm::x");
    let eps = symbol!("canonical_norm::eps");
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s, t],
        &[Atom::var(r) + Atom::var(t)],
        &[vec![vec![Atom::one()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    )?;
    // r=+2 makes r+t=4 on one sheet, but r=-2 is singular. The formal
    // multi-sheet domain deliberately rejects this exceptional specialization.
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::num(4)), (t, Atom::num(2))]),
    };
    assert!(matches!(
        system.pullback(&path, 2),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn a_rational_path_keeps_coordinate_poles_when_the_connection_is_zero() -> Result<()> {
    let s = symbol!("canonical_coordinate::s");
    let x = symbol!("canonical_coordinate::x");
    let eps = symbol!("canonical_coordinate::eps");
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::one()],
        &[vec![vec![Atom::one()]]],
        vec![],
    )?;
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::one() / Atom::var(x))]),
    };
    let pulled = system.pullback(&path, 2)?;
    assert!(
        pulled
            .system
            .matrices
            .iter()
            .flatten()
            .flatten()
            .all(Atom::is_zero)
    );
    let point = BTreeMap::from([(Atom::var(x), Atom::new())]);
    assert!(
        pulled
            .nonzero_conditions
            .iter()
            .any(|a| substitute(a, &point).together().cancel().is_zero())
    );
    Ok(())
}

#[test]
fn canonical_constructor_rejects_nonconstant_matrices_and_invalid_root_letters() {
    let s = symbol!("canonical_input::s");
    let r = symbol!("canonical_input::r");
    let eps = symbol!("canonical_input::eps");
    let root = SquareRoot {
        symbol: r,
        radicand: Atom::var(s),
    };
    assert!(
        CanonicalAlgebraicSystem::new(
            eps,
            &[s],
            &[Atom::var(s)],
            &[vec![vec![Atom::var(eps)]]],
            vec![]
        )
        .is_err()
    );
    assert!(
        CanonicalAlgebraicSystem::new(
            eps,
            &[s],
            &[Atom::var(s)],
            &[vec![vec![Atom::var(s)]]],
            vec![]
        )
        .is_err()
    );
    assert!(
        CanonicalAlgebraicSystem::new(
            eps,
            &[s, s],
            &[Atom::var(s)],
            &[vec![vec![Atom::one()]]],
            vec![]
        )
        .is_err()
    );
    assert!(
        CanonicalAlgebraicSystem::new(
            eps,
            &[s],
            &[Atom::var(r).pow(2) - Atom::var(s)],
            &[vec![vec![Atom::one()]]],
            vec![root.clone()]
        )
        .is_err()
    );
    assert!(
        CanonicalAlgebraicSystem::new(
            eps,
            &[s],
            &[Atom::one() / (Atom::var(r).pow(2) - Atom::var(s))],
            &[vec![vec![Atom::one()]]],
            vec![root]
        )
        .is_err()
    );
}

#[test]
fn removable_denominators_inside_letters_and_radicands_survive_normalization() -> Result<()> {
    let s = symbol!("canonical_removable::s");
    let r = symbol!("canonical_removable::r");
    let x = symbol!("canonical_removable::x");
    let eps = symbol!("canonical_removable::eps");
    let base = Atom::var(s) - 1;
    let expression = Atom::one() / &base + (Atom::var(s).pow(2) - 2) / &base;
    assert!(
        (&expression - Atom::var(s) - 1)
            .together()
            .cancel()
            .is_zero()
    );
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::one())]),
    };
    let letter = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        std::slice::from_ref(&expression),
        &[vec![vec![Atom::one()]]],
        vec![],
    )?;
    assert!(matches!(
        letter.pullback(&path, 1),
        Err(Error::InvalidInput(_))
    ));
    let radicand = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::one()],
        &[vec![vec![Atom::new()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: expression,
        }],
    )?;
    assert!(matches!(
        radicand.pullback(&path, 1),
        Err(Error::InvalidInput(_))
    ));
    let root_base = Atom::var(r) - 1;
    let root_expression = Atom::one() / &root_base + (Atom::var(r).pow(2) - 2) / &root_base;
    let algebraic_letter = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[root_expression],
        &[vec![vec![Atom::one()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    )?;
    assert!(matches!(
        algebraic_letter.pullback(&path, 1),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn removable_path_coordinate_denominators_survive_zero_connections() -> Result<()> {
    let s = symbol!("canonical_removable_path::s");
    let x = symbol!("canonical_removable_path::x");
    let eps = symbol!("canonical_removable_path::eps");
    let base = Atom::var(x) - 1;
    let coordinate = Atom::one() / &base + (Atom::var(x).pow(2) - 2) / &base;
    assert!(
        (&coordinate - Atom::var(x) - 1)
            .together()
            .cancel()
            .is_zero()
    );
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::one()],
        &[vec![vec![Atom::one()]]],
        vec![],
    )?;
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, coordinate)]),
    };
    let pulled = system.pullback(&path, 1)?;
    assert!(
        pulled
            .system
            .matrices
            .iter()
            .flatten()
            .flatten()
            .all(Atom::is_zero)
    );
    let point = BTreeMap::from([(Atom::var(x), Atom::one())]);
    assert!(
        pulled
            .nonzero_conditions
            .iter()
            .any(|a| substitute(a, &point).together().cancel().is_zero())
    );
    Ok(())
}

#[test]
fn complex_root_letters_match_the_analytic_logarithm_series() -> Result<()> {
    use symbolica_amflow::algebraic::RootSeed;
    use symbolica_amflow::diffexp::EpsilonBoundary;
    use symbolica_amflow::{ComplexFloat, FlowOptions, Precision, RunContext};
    let s = symbol!("canonical_complex_letter::s");
    let r = symbol!("canonical_complex_letter::r");
    let x = symbol!("canonical_complex_letter::x");
    let eps = symbol!("canonical_complex_letter::eps");
    let imaginary = Atom::num(Complex::new(Rational::from(0), Rational::from(1)));
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::var(r) + imaginary],
        &[vec![vec![Atom::one()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    )?;
    let pulled = system.pullback(
        &KinematicPath {
            parameter: x,
            coordinates: BTreeMap::from([(s, Atom::var(x))]),
        },
        3,
    )?;
    for (digits, order) in [(50, 48), (70, 80)] {
        let p = Precision::decimal(digits)?;
        let result = pulled.compile(p)?.transport(
            &EpsilonBoundary {
                point: p.i(4),
                leading: 0,
                coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
            },
            &[p.i(9)],
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &FlowOptions {
                digits: 30,
                guard_digits: digits - 30,
                series_order: order,
                ..Default::default()
            },
            &RunContext::default(),
            false,
        )?;
        let g = p.log(&p.div(
            &ComplexFloat::new(p.real(3), p.real(1)),
            &ComplexFloat::new(p.real(2), p.real(1)),
        ));
        let mut expected = p.i(1);
        for (k, row) in result.solution.coefficients.iter().enumerate() {
            if k > 0 {
                expected = p.scale(&p.mul(&expected, &g), 1, k as i64);
            }
            assert!(p.close(&row[0], &expected, 30));
        }
    }
    Ok(())
}
