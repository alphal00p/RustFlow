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
