use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{Error, family::substitute, kinematics::*};

fn equal(a: &Atom, b: &Atom) {
    assert!((a - b).together().cancel().is_zero(), "{a} != {b}");
}

#[test]
fn multivariate_pullback_obeys_the_exact_chain_rule() {
    let (s, t, x, epsilon) = symbol!(
        "path_chain::s",
        "path_chain::t",
        "path_chain::x",
        "path_chain::eps"
    );
    let potential = (Atom::var(s).pow(2) + Atom::var(t)) / (Atom::var(s) - Atom::var(t));
    let regulator = Atom::var(epsilon);
    let system = KinematicSystem {
        epsilon,
        derivatives: [s, t]
            .into_iter()
            .map(|v| (v, vec![vec![&regulator * potential.derivative(v)]]))
            .collect(),
    };
    let parameter = Atom::var(x);
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([
            (s, (Atom::num(1) + &parameter) / (Atom::num(1) - &parameter)),
            (t, Atom::num(2) + parameter.pow(2)),
        ]),
    };
    let pulled = system.pullback(&path).unwrap();
    let restricted_potential = substitute(
        &potential,
        &path
            .coordinates
            .iter()
            .map(|(&v, a)| (Atom::var(v), a.clone()))
            .collect(),
    );
    assert_eq!(pulled.variable, x);
    equal(
        &pulled.matrix[0][0],
        &(&regulator * restricted_potential.derivative(x)),
    );
    assert!(
        pulled.matrix[0][0]
            .get_all_symbols(false)
            .iter()
            .all(|s| [x, epsilon].contains(s))
    );
}

#[test]
fn affine_complex_pullback_preserves_polynomial_numerators_and_all_epsilon_orders() {
    use symbolica::domains::float::Complex;
    use symbolica_amflow::diffexp::EpsilonSystem;

    let (s, t, x, epsilon) = symbol!(
        "affine_polynomial::s",
        "affine_polynomial::t",
        "affine_polynomial::x",
        "affine_polynomial::eps"
    );
    let potential =
        (Atom::var(s) + Atom::num(2) * Atom::var(t)).pow(8) / (Atom::var(s) - Atom::var(t));
    let regulator = Atom::one() + Atom::var(epsilon) + Atom::var(epsilon).pow(7);
    let system = KinematicSystem {
        epsilon,
        derivatives: [s, t]
            .into_iter()
            .map(|v| (v, vec![vec![&regulator * potential.derivative(v)]]))
            .collect(),
    };
    let imaginary = Atom::num(Complex::new(Rational::zero(), Rational::one()));
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, 2 + &imaginary * Atom::var(x)), (t, 1 - Atom::var(x))]),
    };
    let rules = path
        .coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    let expected = substitute(&potential, &rules).derivative(x);
    let pulled = system.pullback(&path).unwrap();
    equal(&pulled.matrix[0][0], &(&regulator * &expected));
    let expanded = EpsilonSystem::from_differential_system(&pulled, epsilon, 7).unwrap();
    let complex_coefficient = BTreeMap::from([(
        Atom::var(symbol!("symbolica_amflow::imaginary_unit")),
        imaginary,
    )]);
    for (k, matrix) in expanded.matrices.iter().enumerate() {
        if [0, 1, 7].contains(&k) {
            equal(&substitute(&matrix[0][0], &complex_coefficient), &expected);
        } else {
            assert!(matrix[0][0].is_zero());
        }
    }
}

#[test]
fn affine_pullback_keeps_weighted_epsilon_pole_cancellation() {
    use symbolica_amflow::diffexp::EpsilonSystem;
    let (s, t, x, epsilon) = symbol!(
        "affine_cancel::s",
        "affine_cancel::t",
        "affine_cancel::x",
        "affine_cancel::eps"
    );
    let eps = Atom::var(epsilon);
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (
                s,
                vec![vec![Atom::one() / (&eps * (1 - &eps * Atom::var(s)))]],
            ),
            (t, vec![vec![-Atom::one() / (2 * &eps)]]),
        ]),
    };
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, 1 + Atom::var(x)), (t, 3 + 2 * Atom::var(x))]),
    };
    let pulled = system.pullback(&path).unwrap();
    let expanded = EpsilonSystem::from_differential_system(&pulled, epsilon, 3).unwrap();
    for (k, matrix) in expanded.matrices.iter().enumerate() {
        equal(
            &matrix[0][0],
            &(Atom::one() + Atom::var(x)).pow(k as i64 + 1),
        );
    }
}

#[test]
fn affine_stationary_coordinates_keep_original_poles_even_at_high_epsilon_order() {
    let (s, t, x, epsilon) = symbol!(
        "affine_stationary::s",
        "affine_stationary::t",
        "affine_stationary::x",
        "affine_stationary::eps"
    );
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::zero()), (t, 1 + Atom::var(x))]),
    };
    for numerator in [Atom::one(), Atom::var(epsilon).pow(8)] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![numerator / Atom::var(s)]]),
                (t, vec![vec![Atom::var(epsilon) / (1 + Atom::var(t))]]),
            ]),
        };
        assert!(matches!(
            system.pullback(&path),
            Err(Error::InvalidInput(_))
        ));
    }
}

#[test]
fn affine_expansion_does_not_truncate_or_accept_nonregular_epsilon_dependence() {
    use symbolica_amflow::diffexp::EpsilonSystem;
    let (s, x, epsilon) = symbol!(
        "affine_epsilon::s",
        "affine_epsilon::x",
        "affine_epsilon::eps"
    );
    let eps = Atom::var(epsilon);
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, 1 + Atom::var(x))]),
    };
    for dependence in [
        eps.pow(-1),
        eps.pow(Atom::num((-1, 2))),
        Atom::one() / (1 - &eps) + eps.pow(Atom::num((101, 2))),
    ] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![dependence / (1 + Atom::var(s))]])]),
        };
        let pulled = system.pullback(&path).unwrap();
        assert!(matches!(
            EpsilonSystem::from_differential_system(&pulled, epsilon, 0),
            Err(Error::Unsupported(_))
        ));
    }
}

#[test]
fn straight_line_retains_exact_complex_endpoints_and_orientation() {
    let (s, x, epsilon) = symbol!("path_line::s", "path_line::x", "path_line::eps");
    let start = BTreeMap::from([(s, Atom::num((1, 3)))]);
    let end = BTreeMap::from([(
        s,
        Atom::num(symbolica::domains::float::Complex::new(
            Rational::from((7, 11)),
            Rational::from((2, 5)),
        )),
    )]);
    let path = KinematicPath::straight_line(x, &start, &end).unwrap();
    for (at, expected) in [(0, &start[&s]), (1, &end[&s])] {
        equal(
            &substitute(
                &path.coordinates[&s],
                &BTreeMap::from([(Atom::var(x), Atom::num(at))]),
            ),
            expected,
        );
    }
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / Atom::var(s)]])]),
    };
    let forward = system.pullback(&path).unwrap();
    let backward = system
        .pullback(&KinematicPath::straight_line(x, &end, &start).unwrap())
        .unwrap();
    let reversed = substitute(
        &forward.matrix[0][0],
        &BTreeMap::from([(Atom::var(x), Atom::num(1) - Atom::var(x))]),
    );
    equal(&backward.matrix[0][0], &-reversed);
}

#[test]
fn canonical_dlog_with_algebraic_letters_matches_direct_logarithmic_derivatives() {
    let (s, t, x, epsilon) = symbol!(
        "path_dlog::s",
        "path_dlog::t",
        "path_dlog::x",
        "path_dlog::eps"
    );
    let letters = [
        Atom::var(s),
        Atom::var(t) - Atom::var(s),
        (Atom::var(s).pow(2) + Atom::var(t)).pow(Atom::num((1, 2))),
    ];
    let matrices = [
        vec![
            vec![Atom::num(1), Atom::num(2)],
            vec![Atom::new(), Atom::num(-1)],
        ],
        vec![
            vec![Atom::new(), Atom::num(3)],
            vec![Atom::num(1), Atom::new()],
        ],
        vec![
            vec![Atom::num((1, 7)), Atom::new()],
            vec![Atom::new(), Atom::num(1)],
        ],
    ];
    let system = KinematicSystem::canonical_dlog(epsilon, &[s, t], &letters, &matrices).unwrap();
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([
            (s, Atom::num(1) + Atom::var(x)),
            (t, Atom::num(3) + Atom::var(x).pow(2)),
        ]),
    };
    let pulled = system.pullback(&path).unwrap();
    let rules = path
        .coordinates
        .iter()
        .map(|(&v, a)| (Atom::var(v), a.clone()))
        .collect();
    for (i, row) in pulled.matrix.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            let expected =
                letters
                    .iter()
                    .zip(&matrices)
                    .fold(Atom::new(), |sum, (letter, coefficients)| {
                        let restricted = substitute(letter, &rules);
                        sum + Atom::var(epsilon) * &coefficients[i][j] * restricted.derivative(x)
                            / restricted
                    });
            equal(entry, &expected);
        }
    }
}

#[test]
fn pullback_rejects_missing_variables_parameter_collisions_and_singular_lines() {
    let (s, t, x, epsilon) = symbol!(
        "path_invalid::s",
        "path_invalid::t",
        "path_invalid::x",
        "path_invalid::eps"
    );
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::var(t).pow(-1)]]),
            (t, vec![vec![Atom::new()]]),
        ]),
    };
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(x)), (t, Atom::num(1))]),
    };
    system.pullback(&path).unwrap();
    let missing = KinematicPath {
        coordinates: BTreeMap::from([(s, Atom::var(x))]),
        ..path.clone()
    };
    assert!(matches!(
        system.pullback(&missing),
        Err(Error::InvalidInput(_))
    ));
    let physical_parameter = KinematicPath {
        parameter: s,
        coordinates: BTreeMap::from([(s, Atom::var(s)), (t, Atom::num(1))]),
    };
    assert!(matches!(
        system.pullback(&physical_parameter),
        Err(Error::InvalidInput(_))
    ));
    let regulator_parameter = KinematicPath {
        parameter: epsilon,
        coordinates: BTreeMap::from([(s, Atom::var(epsilon)), (t, Atom::num(1))]),
    };
    assert!(matches!(
        system.pullback(&regulator_parameter),
        Err(Error::InvalidInput(_))
    ));
    let singular = KinematicPath {
        coordinates: BTreeMap::from([(s, Atom::var(x)), (t, Atom::new())]),
        ..path
    };
    assert!(matches!(
        system.pullback(&singular),
        Err(Error::InvalidInput(_))
    ));
}

#[test]
fn validation_rejects_undeclared_symbols_inexact_numbers_and_function_calls() {
    let (s, t, x, epsilon) = symbol!(
        "path_expr::s",
        "path_expr::t",
        "path_expr::x",
        "path_expr::eps"
    );
    for bad in [
        Atom::var(t),
        Atom::num(Float::with_val(80, 1)),
        parse!("unknown_function(1)"),
    ] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![bad]])]),
        };
        assert!(system.validate().is_err());
    }
    let cross_referenced = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::var(t)), (t, Atom::num(1))]),
    };
    assert!(cross_referenced.validate().is_err());
    let regulator_variable = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(epsilon, vec![vec![Atom::new()]])]),
    };
    assert!(regulator_variable.validate().is_err());
    let dimension_mismatch = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::new()]]),
            (t, vec![vec![Atom::new(); 2]; 2]),
        ]),
    };
    assert!(dimension_mismatch.validate().is_err());
}

#[test]
fn canonical_constructor_and_line_endpoints_validate_declared_data() {
    let (s, t, x, epsilon) = symbol!(
        "path_constructor::s",
        "path_constructor::t",
        "path_constructor::x",
        "path_constructor::eps"
    );
    let letter = [Atom::var(s)];
    let constant = [vec![vec![Atom::num(1)]]];
    assert!(KinematicSystem::canonical_dlog(epsilon, &[s, s], &letter, &constant).is_err());
    assert!(KinematicSystem::canonical_dlog(epsilon, &[epsilon], &letter, &constant).is_err());
    assert!(KinematicSystem::canonical_dlog(epsilon, &[s], &[Atom::new()], &constant).is_err());
    assert!(
        KinematicSystem::canonical_dlog(epsilon, &[s], &letter, &[vec![vec![Atom::var(s)]]])
            .is_err()
    );
    assert!(KinematicSystem::canonical_dlog(epsilon, &[s], &letter, &[]).is_err());
    assert!(
        KinematicPath::straight_line(
            x,
            &BTreeMap::from([(s, Atom::num(1))]),
            &BTreeMap::from([(t, Atom::num(2))])
        )
        .is_err()
    );
    assert!(
        KinematicPath::straight_line(
            x,
            &BTreeMap::from([(s, Atom::var(x))]),
            &BTreeMap::from([(s, Atom::num(2))])
        )
        .is_err()
    );
}
