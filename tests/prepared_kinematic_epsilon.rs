use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{
    diffexp::EpsilonSystem,
    kinematics::{KinematicPath, KinematicSystem, PreparedKinematicEpsilon},
    transport_cache::*,
    *,
};

fn symbols() -> (Symbol, Symbol, Symbol, Symbol) {
    symbol!(
        "prepared_epsilon::s",
        "prepared_epsilon::t",
        "prepared_epsilon::x",
        "prepared_epsilon::eps"
    )
}

fn decode(a: &Atom) -> Atom {
    a.replace(Atom::var(symbol!("symbolica_amflow::imaginary_unit")))
        .with(Atom::i())
}

fn assert_formal_complex_coefficients(a: AtomView<'_>) {
    match a {
        AtomView::Num(n) => match n.get_coeff_view().to_owned() {
            Coefficient::Complex(c) => assert!(c.im.is_zero()),
            _ => panic!("expected exact rational coefficients"),
        },
        AtomView::Var(_) => {}
        AtomView::Add(v) => v.iter().for_each(assert_formal_complex_coefficients),
        AtomView::Mul(v) => v.iter().for_each(assert_formal_complex_coefficients),
        AtomView::Pow(v) => {
            let (base, exponent) = v.get_base_exp();
            assert_formal_complex_coefficients(base);
            assert_formal_complex_coefficients(exponent);
        }
        _ => panic!("expected rational arithmetic"),
    }
}

fn equal_coefficients(legacy: &EpsilonSystem, prepared: &EpsilonSystem) {
    assert_eq!(legacy.variable, prepared.variable);
    assert_eq!(legacy.matrices.len(), prepared.matrices.len());
    for (a, b) in legacy.matrices.iter().zip(&prepared.matrices) {
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().flatten().zip(b.iter().flatten()) {
            // Both APIs use the existing formal imaginary-unit convention.
            assert!(
                (decode(a) - decode(b)).together().cancel().is_zero(),
                "{a} != {b}"
            );
        }
    }
}

#[test]
fn exact_prepared_coefficients_match_legacy_on_rational_and_complex_affine_paths() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps) = (Atom::var(s), Atom::var(t), Atom::var(epsilon));
    let rational: Atom = (1 + &eps + eps.pow(7)) / ((1 - &eps * &sv) * (&sv - &tv));
    let cancelled = ((&sv + &eps).pow(2) - sv.pow(2)) / (&eps * (&sv + &tv));
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (
                s,
                vec![
                    vec![rational.clone(), cancelled],
                    vec![eps.pow(4) / &tv, Atom::i() * &eps / (1 + &eps)],
                ],
            ),
            (
                t,
                vec![
                    vec![-rational / 2, eps.pow(2) / (1 + &eps * &tv)],
                    vec![Atom::zero(), Atom::zero()],
                ],
            ),
        ]),
    };
    let prepared = system.prepare_epsilon(7)?;
    assert_eq!(prepared.max_order(), 7);
    for imaginary in [Atom::zero(), Atom::i()] {
        let path = KinematicPath::straight_line(
            x,
            &BTreeMap::from([(s, Atom::num(3)), (t, Atom::one())]),
            &BTreeMap::from([(s, 4 + &imaginary), (t, 2 - &imaginary)]),
        )?;
        for order in [0, 2, 7] {
            equal_coefficients(
                &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, order)?,
                &prepared.pullback(&path, order)?,
            );
        }
    }
    // Weighted partials cancel exactly, including all retained epsilon orders.
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, 3 + Atom::var(x)), (t, 1 + 2 * Atom::var(x))]),
    };
    assert!(
        prepared
            .pullback(&path, 7)?
            .matrices
            .iter()
            .all(|m| m[0][0].is_zero())
    );
    Ok(())
}

#[test]
fn prepared_native_affine_big_gaussian_rationals_and_weighted_cancellations() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps) = (Atom::var(s), Atom::var(t), Atom::var(epsilon));
    let q =
        parse!("123456789012345678901234567890123456789/987654321098765432109876543210987654319");
    let r =
        parse!("314159265358979323846264338327950288419/271828182845904523536028747135266249777");
    let f: Atom = (&sv + &tv).pow(3) / ((1 + sv.pow(2) + tv.pow(2)) * (1 - &eps * &sv));
    // Along t=2*s+1 this entry becomes t+1. Its removable denominator
    // must cancel in the native output, including at x=1 (s=4).
    let removable = (&tv - 2 * &sv - 1 + (&sv - 4) * (&tv + 1)) / (&sv - 4);
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (
                s,
                vec![
                    vec![2 * &f, (&sv * &tv + Atom::i() * &r) / (&sv + &tv + &eps)],
                    vec![removable * (1 + &eps), Atom::i() * &eps / (&sv + &tv)],
                ],
            ),
            (
                t,
                vec![
                    vec![
                        -&f,
                        (&sv + Atom::i() * &tv) * (1 + &eps + eps.pow(2))
                            / (2 + sv.pow(2) + tv.pow(2)),
                    ],
                    vec![Atom::zero(), -Atom::i() * &eps / (2 * (&sv + &tv))],
                ],
            ),
        ]),
    };
    let prepared = system.prepare_epsilon(4)?;
    // Exercise both source orientations: complex -> physical is the production case.
    let physical = BTreeMap::from([(s, Atom::num(4)), (t, Atom::num(9))]);
    let complex = BTreeMap::from([(s, 3 + Atom::i() * &q), (t, 7 + 2 * Atom::i() * &q)]);
    for (start, end) in [(&complex, &physical), (&physical, &complex)] {
        let path = KinematicPath::straight_line(x, start, end)?;
        for order in [0, 2, 4] {
            let pulled = prepared.pullback(&path, order)?;
            equal_coefficients(
                &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, order)?,
                &pulled,
            );
            for entry in pulled.matrices.iter().flatten().flatten() {
                assert_formal_complex_coefficients(entry.as_view());
            }
            assert!(
                pulled.matrices[0][0][1]
                    .contains_symbol(symbol!("symbolica_amflow::imaginary_unit"))
            );
            // Different Jacobian-weighted partials cancel at every retained order.
            assert!(
                pulled
                    .matrices
                    .iter()
                    .all(|m| m[0][0].is_zero() && m[1][1].is_zero())
            );
            // Native polynomial cancellation must remove the x=1 (or x=0)
            // factor rather than passing a 0/0 expression to the solver.
            let expected = path.coordinates[&s].derivative(x) * (&path.coordinates[&t] + 1);
            assert!(
                (decode(&pulled.matrices[0][1][0]) - expected.expand())
                    .expand()
                    .is_zero()
            );
        }
    }
    Ok(())
}

#[test]
fn prepared_common_denominator_repeated_factors_across_partials_rows_and_orders() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps, xv) = (Atom::var(s), Atom::var(t), Atom::var(epsilon), Atom::var(x));
    let gram: Atom = &sv + &tv;
    let other: Atom = &sv - &tv;
    let weight: Atom = 1 + &eps + eps.pow(2);
    // LCM is 6*gram^3*other^2, rather than a product containing
    // gram^5*other^3. Repeat kernels across both rows and epsilon orders.
    let a: Atom = &weight / (2 * gram.pow(3) * &other);
    let b: Atom = &weight / (3 * gram.pow(2) * other.pow(2));
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![a.clone(), a.clone()], vec![a.clone(), a]]),
            (t, vec![vec![b.clone(), b.clone()], vec![b.clone(), b]]),
        ]),
    };
    let prepared = system.prepare_epsilon(3)?;
    // Two lines through the same prepared snapshot also exercise per-path
    // restriction storage: restrictions from the first line must not survive.
    for slope in [1, 3] {
        let path = KinematicPath {
            parameter: x,
            coordinates: BTreeMap::from([(s, 3 + &xv), (t, 1 + slope * &xv)]),
        };
        for order in [0, 2, 3] {
            let pulled = prepared.pullback(&path, order)?;
            equal_coefficients(
                &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, order)?,
                &pulled,
            );
            let g = &path.coordinates[&s] + &path.coordinates[&t];
            let h = &path.coordinates[&s] - &path.coordinates[&t];
            let expected: Atom = (3 * &h + 2 * slope * &g) / (6 * g.pow(3) * h.pow(2));
            for (k, matrix) in pulled.matrices.iter().enumerate() {
                for entry in matrix.iter().flatten() {
                    if k < 3 {
                        assert!((decode(entry) - &expected).together().cancel().is_zero());
                    } else {
                        assert!(entry.is_zero());
                    }
                }
            }
        }
    }
    Ok(())
}

#[test]
fn prepared_common_denominator_cancels_repeated_factors_only_after_affine_restriction() -> Result<()>
{
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps, xv) = (Atom::var(s), Atom::var(t), Atom::var(epsilon), Atom::var(x));
    let factor: Atom = &sv - 4;
    let denominator: Atom = factor.pow(3);
    let a: Atom = (1 + &eps) / &denominator;
    let b: Atom =
        (1 + &eps) * (-1 + &tv - 2 * &sv - 1 + &denominator * (&tv + 1)) / (2 * &denominator);
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![a]]), (t, vec![vec![b]])]),
    };
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, 1 + &xv), (t, 3 + 2 * &xv)]),
    };
    let pulled = system.prepare_epsilon(2)?.pullback(&path, 2)?;
    equal_coefficients(
        &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 2)?,
        &pulled,
    );
    for k in [0, 1] {
        // Polynomial equality without rational cancellation: all three factors
        // of x-3 must already have disappeared from the returned expression.
        assert!(
            (decode(&pulled.matrices[k][0][0]) - (Atom::num(4) + Atom::num(2) * &xv))
                .expand()
                .is_zero()
        );
        assert!(
            (decode(&pulled.matrices[k][0][0])
                .replace(xv.clone())
                .with(Atom::num(3))
                - 10)
                .is_zero()
        );
    }
    assert!(pulled.matrices[2][0][0].is_zero());
    Ok(())
}

#[test]
fn prepared_common_denominator_preserves_rational_complex_coefficients_and_paths() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps) = (Atom::var(s), Atom::var(t), Atom::var(epsilon));
    let q: Atom =
        parse!("123456789012345678901234567890123456789/987654321098765432109876543210987654319");
    let r: Atom =
        parse!("314159265358979323846264338327950288419/271828182845904523536028747135266249777");
    let shared: Atom = &sv + Atom::i() * &tv + &q;
    let other: Atom = &sv - &tv + &r;
    let c: Atom = &q + Atom::i() * &r;
    let a: Atom = &c * (1 + &eps + eps.pow(2)) / (2 * shared.pow(3) * &other);
    let b: Atom =
        (1 - Atom::i() * &q) * (1 + &eps + eps.pow(2)) / (3 * shared.pow(2) * other.pow(2));
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![a]]), (t, vec![vec![b]])]),
    };
    let prepared = system.prepare_epsilon(3)?;
    let complex = BTreeMap::from([(s, &q + Atom::i() * &r), (t, &r - Atom::i() * &q)]);
    let physical = BTreeMap::from([(s, Atom::num(3)), (t, Atom::one())]);
    for (start, end) in [(&complex, &physical), (&physical, &complex)] {
        let path = KinematicPath::straight_line(x, start, end)?;
        let pulled = prepared.pullback(&path, 3)?;
        equal_coefficients(
            &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 3)?,
            &pulled,
        );
        for entry in pulled.matrices.iter().flatten().flatten() {
            assert_formal_complex_coefficients(entry.as_view());
        }
    }
    Ok(())
}

#[test]
fn prepared_common_denominator_keeps_original_poles_when_all_jacobians_are_zero() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let sv: Atom = Atom::var(s);
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::var(epsilon).pow(8) / sv.pow(3)]]),
            (t, vec![vec![Atom::one() / sv.pow(2)]]),
        ]),
    };
    let prepared = system.prepare_epsilon(0)?;
    for t_coordinate in [Atom::one(), 1 + Atom::var(x)] {
        let path = KinematicPath {
            parameter: x,
            coordinates: BTreeMap::from([(s, Atom::zero()), (t, t_coordinate)]),
        };
        assert!(system.pullback(&path).is_err());
        assert!(matches!(
            prepared.pullback(&path, 0),
            Err(Error::InvalidInput(message)) if message == "path lies in an original kinematic matrix pole"
        ));
    }
    Ok(())
}

// The pre-native-guard check, retained here as an exact rejection oracle for
// explicitly supplied original denominators and epsilon-constant denominators.
fn old_symbolic_prepared_guards(
    path: &KinematicPath,
    original: &[Atom],
    expansion: &[Atom],
) -> Result<()> {
    let substitutions = path
        .coordinates
        .iter()
        .map(|(&s, a)| (Atom::var(s), a.clone()))
        .collect();
    for (poles, expansion) in [(original, false), (expansion, true)] {
        for pole in poles {
            if family::substitute(pole, &substitutions).expand().is_zero() {
                return Err(if expansion {
                    Error::Unsupported("exceptional epsilon specialization".into())
                } else {
                    Error::InvalidInput("original kinematic matrix pole".into())
                });
            }
        }
    }
    Ok(())
}

#[test]
fn prepared_guard_witness_nonzero_source_matches_legacy() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps) = (Atom::var(s), Atom::var(t), Atom::var(epsilon));
    let q =
        parse!("123456789012345678901234567890123456789/987654321098765432109876543210987654319");
    let r =
        parse!("314159265358979323846264338327950288419/271828182845904523536028747135266249777");
    let original = &sv - &tv + &eps * (&sv + &tv);
    let expansion = &sv - &tv;
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![&eps / &original]]),
            (t, vec![vec![Atom::i() * &eps / (2 * &original)]]),
        ]),
    };
    let prepared = system.prepare_epsilon(2)?;
    let complex = BTreeMap::from([(s, &q + Atom::i() * &r), (t, &r - Atom::i() * &q)]);
    let physical = BTreeMap::from([(s, Atom::num(3)), (t, Atom::one())]);
    // At either source, both guards have a nonzero exact epsilon=0 witness.
    for (start, end) in [(&complex, &physical), (&physical, &complex)] {
        let path = KinematicPath::straight_line(x, start, end)?;
        old_symbolic_prepared_guards(
            &path,
            std::slice::from_ref(&original),
            std::slice::from_ref(&expansion),
        )?;
        equal_coefficients(
            &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 2)?,
            &prepared.pullback(&path, 2)?,
        );
    }
    Ok(())
}

#[test]
fn prepared_guard_witness_misses_still_allow_nonconstant_restrictions() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps, xv) = (Atom::var(s), Atom::var(t), Atom::var(epsilon), Atom::var(x));
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, xv.clone()), (t, 1 + &xv)]),
    };
    for original in [
        // The source is an original guard root for every epsilon, but the
        // restriction x is nonconstant. Direct pullback must allow this line.
        sv.clone(),
        // At the source both chosen epsilon values vanish, yet the formal
        // restriction x+epsilon*(epsilon-1)*(1+x) is not identically zero.
        &sv + &eps * (&eps - 1) * &tv,
        // epsilon=0 misses; epsilon=1 alone witnesses the original guard.
        &sv + &eps * &tv,
    ] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![&eps / &original]]),
                (t, vec![vec![Atom::zero()]]),
            ]),
        };
        old_symbolic_prepared_guards(
            &path,
            std::slice::from_ref(&original),
            std::slice::from_ref(&sv),
        )?;
        let pulled = system.prepare_epsilon(3)?.pullback(&path, 3)?;
        equal_coefficients(
            &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 3)?,
            &pulled,
        );
        assert!(!pulled.matrices[1][0][0].is_zero());
    }
    Ok(())
}

#[test]
fn prepared_guard_witness_misses_keep_stationary_and_exceptional_errors() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, tv, eps) = (Atom::var(s), Atom::var(t), Atom::var(epsilon));
    let original = &sv + &eps * (&eps - 1) * &tv;
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            // Above the retained order, and its coordinate is stationary.
            (s, vec![vec![eps.pow(8) / &original]]),
            (t, vec![vec![Atom::zero()]]),
        ]),
    };
    let prepared = system.prepare_epsilon(0)?;
    let exceptional = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::zero()), (t, 1 + Atom::var(x))]),
    };
    // Both witnesses vanish, but the original guard is the nonzero formal
    // polynomial epsilon*(epsilon-1)*(1+x). Only the expansion guard is zero.
    assert!(matches!(
        old_symbolic_prepared_guards(
            &exceptional,
            std::slice::from_ref(&original),
            std::slice::from_ref(&sv),
        ),
        Err(Error::Unsupported(_))
    ));
    assert!(matches!(
        prepared.pullback(&exceptional, 0),
        Err(Error::Unsupported(message))
            if message == "path has an exceptional epsilon specialization; use the legacy pullback"
    ));
    let contained = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::zero()), (t, Atom::zero())]),
    };
    assert!(matches!(
        old_symbolic_prepared_guards(
            &contained,
            std::slice::from_ref(&original),
            std::slice::from_ref(&sv),
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        prepared.pullback(&contained, 0),
        Err(Error::InvalidInput(message))
            if message == "path lies in an original kinematic matrix pole"
    ));
    Ok(())
}

#[test]
fn prepared_native_guards_keep_epsilon_independent_on_large_gaussian_affine_paths() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let (sv, eps, xv) = (Atom::var(s), Atom::var(epsilon), Atom::var(x));
    let q =
        parse!("123456789012345678901234567890123456789/987654321098765432109876543210987654319");
    let r =
        parse!("314159265358979323846264338327950288419/271828182845904523536028747135266249777");
    let shift = &q + Atom::i() * &r;
    let scale = &r - Atom::i() * &q;
    let original = &sv - &shift - &scale * &eps;
    let expansion = &sv - &shift;
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::zero()]]),
            (t, vec![vec![&eps / &original]]),
        ]),
    };
    let prepared = system.prepare_epsilon(3)?;
    let affine = &shift + &scale * &xv;
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, affine.clone()), (t, affine.clone())]),
    };
    old_symbolic_prepared_guards(
        &path,
        std::slice::from_ref(&original),
        std::slice::from_ref(&expansion),
    )?;
    // The restricted original guard is scale*(x-epsilon), not zero. In
    // particular epsilon must not be identified with the path parameter.
    for order in [0, 1, 3] {
        let pulled = prepared.pullback(&path, order)?;
        equal_coefficients(
            &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, order)?,
            &pulled,
        );
        let expected = EpsilonSystem {
            variable: x,
            matrices: (0..=order)
                .map(|k| {
                    vec![vec![if k == 0 {
                        Atom::zero()
                    } else {
                        xv.pow(-(k as i64))
                    }]]
                })
                .collect(),
        };
        equal_coefficients(&expected, &pulled);
    }
    // At s=shift the original guard is -scale*epsilon, which is generically
    // nonzero. The expansion guard vanishes: preserve Unsupported rather than
    // incorrectly classifying this as an original pole by setting epsilon=0.
    // Exercise the same semantics in the checked non-affine fallback too.
    for moving in [affine, &shift + &scale * (&xv + xv.pow(2))] {
        let exceptional_path = KinematicPath {
            parameter: x,
            coordinates: BTreeMap::from([(s, shift.clone()), (t, moving.clone())]),
        };
        assert!(matches!(
            old_symbolic_prepared_guards(
                &exceptional_path,
                std::slice::from_ref(&original),
                std::slice::from_ref(&expansion),
            ),
            Err(Error::Unsupported(_))
        ));
        assert!(matches!(
            prepared.pullback(&exceptional_path, 0),
            Err(Error::Unsupported(message)) if message.contains("exceptional epsilon specialization")
        ));
        // Legacy specializes first and succeeds, giving -dt/dx/scale exactly.
        let legacy = EpsilonSystem::from_differential_system(
            &system.pullback(&exceptional_path)?,
            epsilon,
            0,
        )?;
        equal_coefficients(
            &EpsilonSystem {
                variable: x,
                matrices: vec![vec![vec![-moving.derivative(x) / &scale]]],
            },
            &legacy,
        );
        // A stationary original pole, including formal epsilon and an entry
        // above the retained order, must still reject before its zero Jacobian.
        let stationary_guard = &expansion + &eps * expansion.pow(2);
        let stationary = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![eps.pow(8) / &stationary_guard]]),
                (t, vec![vec![Atom::zero()]]),
            ]),
        };
        assert!(matches!(
            old_symbolic_prepared_guards(
                &exceptional_path,
                std::slice::from_ref(&stationary_guard),
                std::slice::from_ref(&expansion),
            ),
            Err(Error::InvalidInput(_))
        ));
        assert!(stationary.pullback(&exceptional_path).is_err());
        assert!(matches!(
            stationary.prepare_epsilon(0)?.pullback(&exceptional_path, 0),
            Err(Error::InvalidInput(message)) if message.contains("original kinematic matrix pole")
        ));
    }
    Ok(())
}

#[test]
fn prepared_non_affine_fallback_matches_legacy() -> Result<()> {
    let (s, _, x, epsilon) = symbols();
    let system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(
            s,
            vec![vec![(1 + Atom::i() * Atom::var(epsilon)) / Atom::var(s)]],
        )]),
    };
    let prepared = system.prepare_epsilon(2)?;
    let xv = Atom::var(x);
    for coordinate in [
        1 + xv.pow(2),
        Atom::one() / (2 + &xv),
        (Atom::num(2) + &xv).pow(Atom::num((1, 2))),
    ] {
        let path = KinematicPath {
            parameter: x,
            coordinates: BTreeMap::from([(s, coordinate)]),
        };
        equal_coefficients(
            &EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 2)?,
            &prepared.pullback(&path, 2)?,
        );
    }
    Ok(())
}

#[test]
fn prepared_complex_stationary_poles_and_exceptional_epsilon_domains_reject() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let vanishing = Atom::var(s).pow(2) + 1;
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::i()), (t, 2 + Atom::var(x))]),
    };
    let pole = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::var(epsilon).pow(8) / &vanishing]]),
            (t, vec![vec![Atom::zero()]]),
        ]),
    };
    assert!(matches!(
        pole.prepare_epsilon(0)?.pullback(&path, 0),
        Err(Error::InvalidInput(_))
    ));
    let exceptional = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([
            (s, vec![vec![Atom::zero()]]),
            (
                t,
                vec![vec![Atom::var(epsilon) / (&vanishing + Atom::var(epsilon))]],
            ),
        ]),
    };
    assert!(matches!(
        exceptional.prepare_epsilon(0)?.pullback(&path, 0),
        Err(Error::Unsupported(_))
    ));
    let legacy =
        EpsilonSystem::from_differential_system(&exceptional.pullback(&path)?, epsilon, 0)?;
    assert!(legacy.matrices[0][0][0].is_one());
    Ok(())
}

#[test]
fn prepared_depth_and_path_mismatches_are_explicit_and_data_is_a_snapshot() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let mut system = KinematicSystem {
        epsilon,
        derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / Atom::var(s)]])]),
    };
    let prepared = system.prepare_epsilon(2)?;
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, 1 + Atom::var(x))]),
    };
    let before = prepared.pullback(&path, 2)?;
    assert!(matches!(
        prepared.pullback(&path, 3),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(system.prepare_epsilon(1024), Err(Error::Limit(_))));
    for wrong in [
        KinematicPath {
            parameter: x,
            coordinates: BTreeMap::from([(t, Atom::var(x))]),
        },
        KinematicPath {
            parameter: epsilon,
            coordinates: BTreeMap::from([(s, Atom::var(epsilon))]),
        },
        KinematicPath {
            parameter: s,
            coordinates: BTreeMap::from([(s, Atom::var(s))]),
        },
        KinematicPath {
            parameter: symbol!("symbolica_amflow::imaginary_unit"),
            coordinates: BTreeMap::from([(s, Atom::one())]),
        },
    ] {
        assert!(matches!(
            prepared.pullback(&wrong, 1),
            Err(Error::InvalidInput(_))
        ));
    }
    system.derivatives.get_mut(&s).unwrap()[0][0] = Atom::zero();
    equal_coefficients(&before, &prepared.pullback(&path, 2)?);
    fn assert_send_sync<T: Send + Sync>() {}
    assert_send_sync::<KinematicSystem>();
    assert_send_sync::<PreparedKinematicEpsilon>();
    assert_send_sync::<RustFlow>();
    Ok(())
}

#[test]
fn nonregular_partials_and_fractional_epsilon_powers_are_not_truncated_away() -> Result<()> {
    let (s, t, x, epsilon) = symbols();
    let eps = Atom::var(epsilon);
    for entry in [
        eps.pow(-1),
        eps.pow(Atom::num((-1, 2))),
        Atom::one() / (1 - &eps) + eps.pow(Atom::num((101, 2))),
    ] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![entry]])]),
        };
        assert!(matches!(
            system.prepare_epsilon(0),
            Err(Error::Unsupported(_))
        ));
    }
    // Some legacy lines cancel epsilon poles only after weighting partials.
    // Keep that API working; explicit preparation requires regular partials.
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
    assert!(matches!(
        system.prepare_epsilon(1),
        Err(Error::Unsupported(_))
    ));
    let legacy = EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 1)?;
    assert!(
        (legacy.matrices[0][0][0].clone() - (Atom::one() + Atom::var(x)))
            .together()
            .cancel()
            .is_zero()
    );
    Ok(())
}

#[test]
fn exceptional_specializations_and_stationary_original_poles_reject_before_truncation() -> Result<()>
{
    let (s, t, x, epsilon) = symbols();
    let path = KinematicPath {
        parameter: x,
        coordinates: BTreeMap::from([(s, Atom::zero()), (t, 1 + Atom::var(x))]),
    };
    for entry in [
        Atom::one() / Atom::var(s),
        Atom::var(epsilon).pow(8) / Atom::var(s),
    ] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![entry]]), (t, vec![vec![Atom::zero()]])]),
        };
        assert!(system.pullback(&path).is_err());
        assert!(matches!(
            system.prepare_epsilon(0)?.pullback(&path, 0),
            Err(Error::InvalidInput(_))
        ));
    }
    for numerator in [Atom::one(), Atom::var(epsilon)] {
        let system = KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([
                (s, vec![vec![Atom::zero()]]),
                (
                    t,
                    vec![vec![numerator / (Atom::var(s) + Atom::var(epsilon))]],
                ),
            ]),
        };
        assert!(matches!(
            system.prepare_epsilon(0)?.pullback(&path, 0),
            Err(Error::Unsupported(_))
        ));
        // epsilon/(s+epsilon) becomes 1 on this line; the prepared route must
        // reject even this regular legacy specialization instead of returning 0.
        if let Ok(legacy) =
            EpsilonSystem::from_differential_system(&system.pullback(&path)?, epsilon, 0)
        {
            assert!(legacy.matrices[0][0][0].is_one());
        }
    }
    Ok(())
}

fn flow(conditions: &[Atom]) -> Result<RustFlow> {
    let (s, _, _, epsilon) = symbols();
    RustFlow::with_conditions(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / Atom::var(s)]])]),
        },
        &[parse!("prepared_epsilon::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "positive s, principal logarithm",
        conditions,
    )
}

fn seed(flow: &RustFlow, bits: u32, verified: u32) -> Result<RustFlowCache> {
    let p = Precision { bits };
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(symbols().0, Atom::one())])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(-1, 1)?,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            verified,
            bits,
            vec![vec![p.real(0)]; 3],
            "analytic seed, finite declared source evidence",
        )?,
    })?;
    Ok(cache)
}

fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}

#[test]
fn prepared_transport_matches_legacy_at_80_and_201_input_bits_and_keeps_source_floor() -> Result<()>
{
    let legacy = flow(&[])?;
    let prepared = legacy.prepare_epsilon(2)?;
    assert_eq!(legacy.identity().key(), prepared.identity().key());
    let destination = BTreeMap::from([(symbols().0, Atom::num(2))]);
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let options = FlowOptions {
            digits,
            guard_digits: 20,
            series_order: 32,
            ..Default::default()
        };
        let mut old_cache = seed(&legacy, bits, verified)?;
        let mut new_cache = old_cache.clone();
        let range = EpsilonRange::new(-1, 1)?;
        let old = legacy.evaluate_to(
            &mut old_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy(),
        )?;
        let new = prepared.evaluate_to(
            &mut new_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy(),
        )?;
        assert_eq!(old.boundary.coefficients, new.boundary.coefficients);
        assert_eq!(
            old.boundary.accuracy.comparison_errors(),
            new.boundary.accuracy.comparison_errors()
        );
        assert_eq!(new.boundary.accuracy.input_verified_digits(), verified);
        assert!(new.boundary.accuracy.verified_digits() <= verified);
        let p = Precision {
            bits: new.boundary.accuracy.working_bits(),
        };
        assert!(p.close(&new.boundary.coefficients[1][0], &p.log(&p.i(2)), digits));
        assert!(
            new.boundary
                .accuracy
                .comparison_errors()
                .iter()
                .flatten()
                .all(|e| *e >= p.tolerance(verified))
        );
        let hit = prepared.evaluate_to(
            &mut new_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy(),
        )?;
        assert!(hit.transport.is_none());
        let shallow = legacy.prepare_epsilon(0)?;
        let size = new_cache.len();
        assert!(matches!(
            shallow.evaluate_to(
                &mut new_cache,
                &destination,
                range,
                &options,
                &RunContext::default(),
                &policy()
            ),
            Err(Error::InvalidInput(_))
        ));
        assert_eq!(new_cache.len(), size);
    }
    Ok(())
}

#[test]
fn prepared_complex_source_to_physical_transport_matches_legacy_at_80_and_201_bits() -> Result<()> {
    let legacy = flow(&[])?;
    let prepared = legacy.prepare_epsilon(2)?;
    let source: Atom = 1 + Atom::i() * Atom::num((1, 100));
    let destination = BTreeMap::from([(symbols().0, Atom::num(2))]);
    for (bits, verified, digits) in [(80, 16, 10), (201, 48, 32)] {
        let options = FlowOptions {
            digits,
            guard_digits: 20,
            series_order: 32,
            ..Default::default()
        };
        let mut old_cache = RustFlowCache::default();
        let boundary = seed(&legacy, bits, verified)?
            .entries()
            .first()
            .unwrap()
            .clone();
        old_cache.insert(CachedBoundary {
            point: CachedPoint::Exact(BTreeMap::from([(symbols().0, source.clone())])),
            ..boundary
        })?;
        let mut new_cache = old_cache.clone();
        let range = EpsilonRange::new(-1, 1)?;
        let old = legacy.evaluate_to(
            &mut old_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy(),
        )?;
        let new = prepared.evaluate_to(
            &mut new_cache,
            &destination,
            range,
            &options,
            &RunContext::default(),
            &policy(),
        )?;
        let p = Precision {
            bits: new.boundary.accuracy.working_bits(),
        };
        for (a, b) in old
            .boundary
            .coefficients
            .iter()
            .flatten()
            .zip(new.boundary.coefficients.iter().flatten())
        {
            assert!(p.close(a, b, digits));
        }
        assert_eq!(new.boundary.accuracy.input_verified_digits(), verified);
        assert!(new.boundary.accuracy.verified_digits() <= verified);
        assert!(
            new.boundary
                .accuracy
                .comparison_errors()
                .iter()
                .flatten()
                .all(|e| *e >= p.tolerance(verified))
        );
        // Y(source)=1/epsilon, so Y(target)=(target/source)^epsilon/epsilon.
        let expected = p.log(&p.div(&p.i(200), &p.complex(100, 1)));
        assert!(p.close(&new.boundary.coefficients[1][0], &expected, digits));
        assert!(p.close(
            &new.boundary.coefficients[2][0],
            &p.div(&p.mul(&expected, &expected), &p.i(2)),
            digits
        ));
    }
    Ok(())
}

#[test]
fn prepared_flow_preserves_domain_precision_and_cancellation_rejections() -> Result<()> {
    let guarded = flow(&[Atom::var(symbols().0) - 2])?;
    let prepared = guarded.prepare_epsilon(2)?;
    let range = EpsilonRange::new(-1, 1)?;
    let options = FlowOptions {
        digits: 10,
        guard_digits: 20,
        ..Default::default()
    };
    let mut cache = seed(&guarded, 80, 16)?;
    for destination in [0, 2, 3] {
        assert!(
            prepared
                .evaluate_to(
                    &mut cache,
                    &BTreeMap::from([(symbols().0, Atom::num(destination))]),
                    range,
                    &options,
                    &RunContext::default(),
                    &policy()
                )
                .is_err()
        );
        assert_eq!(cache.len(), 1);
    }
    let mut weak = seed(&guarded, 201, 8)?;
    assert!(
        prepared
            .evaluate_to(
                &mut weak,
                &BTreeMap::from([(symbols().0, Atom::num((3, 2)))]),
                range,
                &options,
                &RunContext::default(),
                &policy()
            )
            .is_err()
    );
    assert_eq!(weak.len(), 1);
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        prepared.evaluate_to(
            &mut cache,
            &BTreeMap::from([(symbols().0, Atom::num((3, 2)))]),
            range,
            &options,
            &context,
            &policy()
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(cache.len(), 1);
    let unguarded = flow(&[])?.prepare_epsilon(2)?;
    let mut unrelated = seed(&unguarded, 80, 16)?;
    assert!(matches!(
        prepared.evaluate_to(
            &mut unrelated,
            &BTreeMap::from([(symbols().0, Atom::one())]),
            range,
            &options,
            &RunContext::default(),
            &policy()
        ),
        Err(Error::IncompleteReduction(_))
    ));
    assert_eq!(unrelated.len(), 1);
    Ok(())
}

#[test]
fn prepared_flow_still_checks_diagonal_gauge_against_the_original_connection() -> Result<()> {
    let legacy = flow(&[])?;
    let prepared = legacy.prepare_epsilon(2)?;
    let options = FlowOptions {
        digits: 10,
        guard_digits: 20,
        diagonal_gauge_estimate: Some(physical_transport::DiagonalGaugeEstimate {
            factors: vec![Atom::var(symbols().0)],
            radicands: vec![Atom::one()],
        }),
        ..Default::default()
    };
    for flow in [&legacy, &prepared] {
        let mut cache = seed(flow, 80, 16)?;
        assert!(
            matches!(flow.evaluate_to(&mut cache, &BTreeMap::from([(symbols().0, Atom::num(2))]), EpsilonRange::new(-1, 1)?, &options, &RunContext::default(), &policy()), Err(Error::InvalidInput(message)) if message.contains("logarithmic derivative"))
        );
        assert_eq!(cache.len(), 1);
    }
    Ok(())
}
