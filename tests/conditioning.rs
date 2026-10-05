use std::cell::RefCell;
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicSystem, RootSeed, SquareRoot};
use symbolica_amflow::diffexp::{EpsilonBoundary, EpsilonSystem, transport_epsilon};
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 40,
        series_order: 16,
        ..Default::default()
    }
}
fn source() -> EpsilonSystem {
    EpsilonSystem {
        variable: symbol!("conditioning::x"),
        matrices: vec![vec![
            vec![Atom::zero(), Atom::one() - Atom::num(10).pow(100)],
            vec![Atom::zero(), Atom::zero()],
        ]],
    }
}
fn boundary(p: Precision) -> Result<EpsilonBoundary> {
    Ok(EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![
            p.eval(&Atom::num(10).pow(100), &Default::default())?,
            p.i(1),
        ]],
    })
}

#[test]
fn independent_zero_profiles_trigger_a_precision_jump_then_two_successes() -> Result<()> {
    let calls = RefCell::new(Vec::new());
    let result = transport_epsilon(
        &source(),
        |p| {
            calls.borrow_mut().push(p.bits);
            boundary(p)
        },
        &[Atom::one()],
        &options(),
        &RunContext::default(),
        true,
    )?;
    let p = Precision {
        bits: result.diagnostics.working_bits,
    };
    assert!(p.close(&result.coefficients[0][0], &p.i(1), 20));
    assert_eq!(result.verified_digits, Some(20));
    let calls = calls.into_inner();
    assert_eq!(calls.len(), 3, "{calls:?}");
    assert!(calls[1] > Precision::decimal(120)?.bits);
    assert!(calls[2] > calls[1]);
    let checked = result.diagnostics.conditioning_digits.unwrap();
    assert!(
        checked >= options().digits && checked < p.bits / 4,
        "{checked}"
    );
    assert!(result.checkpoints.iter().any(|c| c.point == p.i(1)));
    Ok(())
}
#[test]
fn fixed_precision_and_exhausted_refinement_are_typed_failures() -> Result<()> {
    let p = Precision::decimal(60)?;
    let error = source()
        .compile(p, &Default::default())?
        .transport(
            &boundary(p)?,
            &[p.i(1)],
            &options(),
            &RunContext::default(),
            false,
        )
        .unwrap_err();
    assert!(matches!(error,Error::InsufficientPrecision{minimum_bits,..} if minimum_bits>p.bits));
    let mut o = options();
    o.max_precision_attempts = 1;
    let error = transport_epsilon(
        &source(),
        boundary,
        &[Atom::one()],
        &o,
        &RunContext::default(),
        false,
    )
    .unwrap_err();
    assert!(matches!(error, Error::Accuracy(_)));
    Ok(())
}
#[test]
fn exact_waypoint_collapse_retries_but_numeric_zero_length_is_identity() -> Result<()> {
    let x = symbol!("conditioning_coordinates::x");
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![vec![
            vec![Atom::zero(), Atom::one()],
            vec![Atom::zero(), Atom::zero()],
        ]],
    };
    let huge = Atom::num(10).pow(100);
    let result = transport_epsilon(
        &system,
        |p| {
            Ok(EpsilonBoundary {
                point: p.eval(&huge, &Default::default())?,
                leading: 0,
                coefficients: vec![vec![p.zero(), p.i(1)]],
            })
        },
        &[&huge + Atom::one()],
        &options(),
        &RunContext::default(),
        true,
    )?;
    let p = Precision {
        bits: result.diagnostics.working_bits,
    };
    assert!(p.close(&result.coefficients[0][0], &p.i(1), 20));
    assert!(p.bits > Precision::decimal(120)?.bits);
    let p = Precision::decimal(60)?;
    let point = p.eval(&huge, &Default::default())?;
    let identity = system.compile(p, &Default::default())?.transport(
        &EpsilonBoundary {
            point: point.clone(),
            leading: 0,
            coefficients: vec![vec![p.i(7), p.i(1)]],
        },
        std::slice::from_ref(&point),
        &options(),
        &RunContext::default(),
        false,
    )?;
    assert_eq!(identity.coefficients[0][0], p.i(7));
    assert_eq!(identity.diagnostics.steps, 0);
    let exact_numeric = Atom::num(Complex::new(point.re.to_rational(), point.im.to_rational()));
    let identity = transport_epsilon(
        &system,
        |q| {
            Ok(EpsilonBoundary {
                point: point.clone(),
                leading: 0,
                coefficients: vec![vec![q.i(7), q.i(1)]],
            })
        },
        &[exact_numeric],
        &options(),
        &RunContext::default(),
        false,
    )?;
    assert_eq!(identity.diagnostics.steps, 0);
    assert_eq!(identity.coefficients[0][0], p.i(7));

    Ok(())
}

#[test]
fn higher_precision_physical_endpoint_matches_computed_and_saved_values() -> Result<()> {
    let p = Precision::decimal(60)?;
    let coordinates = Precision::decimal(160)?;
    let center = coordinates.powi(&coordinates.i(10), 100);
    let target = coordinates.add(&center, &coordinates.i(1));
    assert_ne!(target, p.add(&center, &p.sub(&target, &center)));
    let system = EpsilonSystem {
        variable: symbol!("conditioning_declared_endpoint::x"),
        matrices: vec![vec![
            vec![Atom::zero(), Atom::one()],
            vec![Atom::zero(), Atom::zero()],
        ]],
    };
    let result = system.compile(p, &Default::default())?.transport(
        &EpsilonBoundary {
            point: center.clone(),
            leading: 0,
            coefficients: vec![vec![p.zero(), p.i(1)]],
        },
        std::slice::from_ref(&target),
        &options(),
        &RunContext::default(),
        true,
    )?;
    assert_eq!(result.point, target);
    assert_eq!(result.coefficients[0][0], p.i(1));
    assert_eq!(result.segments.len(), 1);
    assert_eq!(result.segments[0].end, target);
    assert_eq!(result.evaluate_path(&target)?[0][0], p.i(1));
    let midpoint = coordinates.add(&center, &coordinates.scale(&coordinates.i(1), 1, 2));
    assert_eq!(
        result.evaluate_path(&midpoint)?[0][0],
        p.scale(&p.i(1), 1, 2)
    );
    Ok(())
}
#[test]
fn complex_mapped_root_and_epsilon_representations_cannot_hide_cancellation() -> Result<()> {
    let p = Precision::decimal(60)?;
    let x = symbol!("conditioning_root::x");
    let r = symbol!("conditioning_root::r");
    let huge = Atom::num(10).pow(100);
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![
                vec![
                    Atom::zero(),
                    (Atom::num(Complex::new(Rational::one(), Rational::from(2))) - &huge)
                        * Atom::var(r),
                ],
                vec![Atom::zero(), Atom::zero()],
            ],
        },
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::one(),
        }],
    );
    let error = system
        .compile(p)?
        .transport_ordinary(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.eval(&huge, &Default::default())?, p.i(1)],
            },
            &[p.i(1)],
            &BTreeMap::from([(r, RootSeed::Principal)]),
            &options(),
            &RunContext::default(),
            false,
        )
        .unwrap_err();
    assert!(
        matches!(error, Error::InsufficientPrecision { .. }),
        "{error:?}"
    );
    let high = Precision::decimal(160)?;
    let resolved = system.compile(high)?.transport_ordinary(
        &BoundaryData {
            point: high.zero(),
            values: vec![high.eval(&huge, &Default::default())?, high.i(1)],
        },
        &[high.i(1)],
        &BTreeMap::from([(r, RootSeed::Principal)]),
        &options(),
        &RunContext::default(),
        false,
    )?;
    assert!(high.close(
        &resolved.solution.coefficients[0][0],
        &high.complex(1, 2),
        20
    ));
    // Actual Mobius map, with a decoupled pole supplying the chart geometry.
    let mut system = source();
    system.matrices[0][0].push(Atom::zero());
    system.matrices[0][1].push(Atom::zero());
    system.matrices[0].push(vec![
        Atom::zero(),
        Atom::zero(),
        Atom::one() / (Atom::var(system.variable).pow(2) - Atom::num(9)),
    ]);
    let mut o = options();
    o.local_coordinate = LocalCoordinate::BalancedMobius;
    o.series_order = 64;
    let result = transport_epsilon(
        &system,
        |p| {
            let mut b = boundary(p)?;
            b.coefficients[0].push(p.zero());
            Ok(b)
        },
        &[Atom::one()],
        &o,
        &RunContext::default(),
        true,
    )?;
    let p = Precision {
        bits: result.diagnostics.working_bits,
    };
    assert!(p.close(&result.coefficients[0][0], &p.i(1), 20));
    assert!(result.segments.iter().any(|s| matches!(
        s.coordinate,
        symbolica_amflow::local_coordinates::TaylorCoordinate::Mobius(_)
    )));
    Ok(())
}
#[test]
fn strong_cache_evidence_never_exceeds_the_conditioning_ceiling() -> Result<()> {
    let s = symbol!("conditioning_cache::s");
    let eps = symbol!("conditioning_cache::eps");
    let system = source();
    let flow = RustFlow::new(
        KinematicSystem {
            epsilon: eps,
            derivatives: BTreeMap::from([(s, system.matrices[0].clone())]),
        },
        &[
            parse!("conditioning_cache::I0"),
            parse!("conditioning_cache::I1"),
        ],
        &Atom::one(),
        Prescription::PlusI0,
        "regular exact nilpotent system",
    )?;
    let p = Precision::decimal(250)?;
    let mut cache = RustFlowCache::default();
    cache.insert(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(BTreeMap::from([(s, Atom::zero())])),
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 0)?,
        coefficients: boundary(p)?.coefficients,
        accuracy: BoundaryAccuracy::supplied(
            180,
            p.bits,
            vec![vec![p.real(0); 2]],
            "independent exact integer seed",
        )?,
    })?;
    let result = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::one())]),
        EpsilonRange::new(0, 0)?,
        &options(),
        &RunContext::default(),
        &ScaledDistance {
            scales: BTreeMap::new(),
            admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
        },
    )?;
    let ceiling = result
        .transport
        .as_ref()
        .unwrap()
        .diagnostics
        .conditioning_digits
        .unwrap();
    assert!(result.boundary.accuracy.verified_digits() <= ceiling.saturating_sub(2).max(20));
    assert!(result.boundary.accuracy.verified_digits() < 50);
    assert!(p.close(&result.boundary.coefficients[0][0], &p.i(1), 20));
    assert!(
        cache
            .entries()
            .iter()
            .skip(1)
            .all(|b| b.accuracy.verified_digits() <= ceiling)
    );
    let mut weak = cache.entries()[0].clone();
    weak.accuracy = BoundaryAccuracy::supplied(
        40,
        p.bits,
        vec![vec![p.real(0); 2]],
        "limited input evidence",
    )?;
    let mut weak_bank = RustFlowCache::default();
    weak_bank.insert(weak)?;
    let failure = flow.evaluate_to(
        &mut weak_bank,
        &BTreeMap::from([(s, Atom::one())]),
        EpsilonRange::new(0, 0)?,
        &options(),
        &RunContext::default(),
        &ScaledDistance {
            scales: BTreeMap::new(),
            admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
        },
    );
    assert!(matches!(failure, Err(Error::Accuracy(_))));
    assert_eq!(weak_bank.len(), 1);
    Ok(())
}
#[test]
fn cancellation_and_step_limits_are_not_precision_retries() -> Result<()> {
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        transport_epsilon(
            &source(),
            boundary,
            &[Atom::one()],
            &options(),
            &context,
            false
        ),
        Err(Error::Cancelled)
    ));
    let mut o = options();
    o.max_steps = 1;
    o.series_order = 8;
    let x = symbol!("conditioning_limit::x");
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![vec![vec![Atom::one()]]],
    };
    assert!(matches!(
        transport_epsilon(
            &system,
            |p| Ok(EpsilonBoundary {
                point: p.zero(),
                leading: 0,
                coefficients: vec![vec![p.i(1)]]
            }),
            &[Atom::num(10)],
            &o,
            &RunContext::default(),
            false
        ),
        Err(Error::Limit(_))
    ));
    Ok(())
}

#[test]
fn epsilon_hierarchy_uses_the_same_precision_retry() -> Result<()> {
    let mut system = source();
    let shift = system.matrices.remove(0);
    system.matrices = vec![vec![vec![Atom::zero(); 2]; 2], shift];
    let result = transport_epsilon(
        &system,
        |p| {
            Ok(EpsilonBoundary {
                point: p.zero(),
                leading: 0,
                coefficients: vec![
                    vec![p.zero(), p.i(1)],
                    vec![
                        p.eval(&Atom::num(10).pow(100), &Default::default())?,
                        p.zero(),
                    ],
                ],
            })
        },
        &[Atom::one()],
        &options(),
        &RunContext::default(),
        false,
    )?;
    let p = Precision {
        bits: result.diagnostics.working_bits,
    };
    assert!(p.close(&result.coefficients[1][0], &p.i(1), 20));
    assert_eq!(result.coefficients[0], vec![p.zero(), p.i(1)]);
    Ok(())
}

#[test]
fn requested_accuracy_is_not_reduced_by_a_coarse_ceiling_conversion() -> Result<()> {
    let x = symbol!("conditioning_high_digits::x");
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![vec![vec![Atom::zero()]]],
    };
    let mut o = options();
    o.digits = 100;
    o.guard_digits = 0;
    let result = transport_epsilon(
        &system,
        |p| {
            Ok(EpsilonBoundary {
                point: p.zero(),
                leading: 0,
                coefficients: vec![vec![p.i(1)]],
            })
        },
        &[Atom::one()],
        &o,
        &RunContext::default(),
        false,
    )?;
    assert_eq!(result.verified_digits, Some(100));
    assert!(result.diagnostics.conditioning_digits.unwrap() >= 100);
    Ok(())
}

#[test]
fn saved_path_does_not_hide_an_unresolved_overlapping_segment() -> Result<()> {
    let system = EpsilonSystem {
        variable: symbol!("conditioning_overlap::x"),
        matrices: vec![vec![vec![Atom::zero()]]],
    };
    let mut solution = transport_epsilon(
        &system,
        |p| {
            Ok(EpsilonBoundary {
                point: p.zero(),
                leading: 0,
                coefficients: vec![vec![p.i(1)]],
            })
        },
        &[Atom::one()],
        &options(),
        &RunContext::default(),
        true,
    )?;
    let p = Precision {
        bits: solution.diagnostics.working_bits,
    };
    assert_eq!(solution.evaluate_path(&p.i(1))?[0][0], p.i(1));
    let mut unresolved = solution.segments[0].clone();
    let huge = p.eval(&Atom::num(10).pow(100), &Default::default())?;
    unresolved.coefficients = vec![vec![huge.clone()], vec![p.sub(&p.i(1), &huge)]];
    solution.segments.push(unresolved);
    assert!(matches!(
        solution.evaluate_path(&p.i(1)),
        Err(Error::InsufficientPrecision { .. })
    ));
    Ok(())
}
