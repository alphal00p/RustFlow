use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, RootSeed, SquareRoot};
use symbolica_amflow::diffexp::EpsilonBoundary;
use symbolica_amflow::kinematics::KinematicPath;
use symbolica_amflow::{FlowOptions, Precision, RunContext};

#[test]
fn algebraic_total_derivative_and_path_pullback_match_exactly() {
    let a = symbol!("algebraic_path::a");
    let b = symbol!("algebraic_path::b");
    let r = symbol!("algebraic_path::r");
    let x = symbol!("algebraic_path::x");
    let eps = symbol!("algebraic_path::eps");
    let letter = (Atom::var(a) + Atom::var(r)) / (Atom::var(a) - Atom::var(r));
    let system = AlgebraicKinematicSystem::canonical_dlog(
        eps,
        &[a, b],
        &[letter],
        &[vec![vec![Atom::num(1)]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(b).pow(2) - Atom::var(a),
        }],
    )
    .unwrap();
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([
            (a, Atom::num(1) + Atom::var(x)),
            (b, Atom::num(4) + Atom::var(x)),
        ]),
    };
    let pulled = system.pullback(&path, 3).unwrap();
    assert_eq!(
        pulled.roots[0].radicand,
        (Atom::num(15) + Atom::var(x) * 7 + Atom::var(x).pow(2)).expand()
    );
    let expected =
        (Atom::num(23) + Atom::var(x) * 5) / (Atom::var(r) * (Atom::num(14) + Atom::var(x) * 5));
    // Compare modulo r²=R with the same native quotient reduction used at compile.
    let difference = (&pulled.system.matrices[1][0][0] - expected)
        .together()
        .cancel();
    let reduced = difference
        .replace(Atom::var(r).pow(2))
        .with(pulled.roots[0].radicand.clone());
    assert!(reduced.together().cancel().is_zero());
    assert!(pulled.system.matrices[0][0][0].is_zero());
    let p = Precision::decimal(50).unwrap();
    let compiled = pulled.compile(p).unwrap();
    let boundary = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()], vec![p.zero()]],
    };
    let result = compiled
        .transport(
            &boundary,
            &[p.i(1)],
            &BTreeMap::from([(r, RootSeed::Value(p.i(10)))]),
            &FlowOptions {
                digits: 20,
                guard_digits: 30,
                series_order: 64,
                ..Default::default()
            },
            &RunContext::default(),
            false,
        )
        .unwrap();
    let half = p.rational(&Rational::from((1, 2)));
    let first = p.pow(&p.i(15), &half);
    let last = p.pow(&p.i(23), &half);
    let first_letter = p.div(&p.add(&p.i(1), &first), &p.sub(&p.i(1), &first));
    let last_letter = p.div(&p.add(&p.i(2), &last), &p.sub(&p.i(2), &last));
    let logarithm = p.log(&p.div(&last_letter, &first_letter));
    let mut expected = p.i(1);
    for k in 0..4 {
        assert!(p.norm(&p.sub(&result.solution.coefficients[k][0], &expected)) < p.tolerance(32));
        expected = p.scale(&p.mul(&expected, &logarithm), 1, k as i64 + 1);
    }
}

#[test]
fn algebraic_path_rejects_undeclared_coordinates_and_root_dependent_radicands() {
    let a = symbol!("algebraic_path_errors::a");
    let r = symbol!("algebraic_path_errors::r");
    let x = symbol!("algebraic_path_errors::x");
    let eps = symbol!("algebraic_path_errors::eps");
    let mut system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([(a, vec![vec![Atom::var(eps) / Atom::var(r)]])]),
        roots: vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(a),
        }],
    };
    let wrong = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(r, Atom::var(x))]),
    };
    assert!(system.pullback(&wrong, 2).is_err());
    let colliding = KinematicPath {
        parameter: r,
        coordinates: BTreeMap::from([(a, Atom::var(r))]),
    };
    assert!(system.pullback(&colliding, 2).is_err());
    system.roots[0].radicand = Atom::var(r);
    assert!(system.validate().is_err());
    system.roots[0].radicand = Atom::var(eps);
    assert!(system.validate().is_err());
    assert!(
        AlgebraicKinematicSystem::canonical_dlog(
            eps,
            &[a],
            &[Atom::var(r).pow(2) - Atom::var(a)],
            &[vec![vec![Atom::num(1)]]],
            vec![SquareRoot {
                symbol: r,
                radicand: Atom::var(a)
            }]
        )
        .is_err()
    );
}

#[test]
fn stationary_coordinates_cannot_erase_source_or_root_poles() {
    let s = symbol!("algebraic_path_poles::s");
    let t = symbol!("algebraic_path_poles::t");
    let r = symbol!("algebraic_path_poles::r");
    let x = symbol!("algebraic_path_poles::x");
    let eps = symbol!("algebraic_path_poles::eps");
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::new()), (t, Atom::var(x))]),
    };
    let mut system = AlgebraicKinematicSystem {
        epsilon: eps,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::num(1) / Atom::var(s)]]),
            (t, vec![vec![Atom::new()]]),
        ]),
        roots: vec![],
    };
    assert!(matches!(
        system.pullback(&path, 2),
        Err(symbolica_amflow::Error::InvalidInput(_))
    ));
    system.derivatives.get_mut(&s).unwrap()[0][0] = Atom::new();
    system.roots = vec![SquareRoot {
        symbol: r,
        radicand: Atom::num(1) / Atom::var(s),
    }];
    assert!(matches!(
        system.pullback(&path, 2),
        Err(symbolica_amflow::Error::InvalidInput(_))
    ));
    system.roots[0].radicand = Atom::var(s);
    assert!(matches!(
        system.pullback(&path, 2),
        Err(symbolica_amflow::Error::InvalidInput(_))
    ));
    system.roots[0].radicand = Atom::var(t) + 1;
    system.derivatives.get_mut(&s).unwrap()[0][0] =
        Atom::num(1) / (Atom::var(r).pow(2) - Atom::var(t) - 1 + Atom::var(s));
    assert!(matches!(
        system.pullback(&path, 2),
        Err(symbolica_amflow::Error::InvalidInput(_))
    ));
    system.derivatives.get_mut(&s).unwrap()[0][0] = Atom::num(1) / (Atom::var(t) + 1);
    assert!(system.pullback(&path, 2).is_ok());
}
