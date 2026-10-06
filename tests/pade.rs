use symbolica::prelude::*;
use symbolica_amflow::{
    diffexp::{EpsilonBoundary, EpsilonSystem, transport_epsilon},
    local_coordinates::LocalCoordinate,
    *,
};

fn options(degree: usize, order: usize, working: u32) -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: working - 20,
        series_order: order,
        pade: Some(PadeOptions {
            degree,
            ..Default::default()
        }),
        max_steps: 10000,
        ..Default::default()
    }
}
fn rational() -> EpsilonSystem {
    let x = symbol!("pade_integration::x");
    EpsilonSystem {
        variable: x,
        matrices: vec![vec![vec![-Atom::one() / (Atom::one() + Atom::var(x))]]],
    }
}
#[test]
fn accepted_rational_values_and_saved_samples_share_representation() -> Result<()> {
    // A degree-one rational is exact at low order, whereas its order-eight
    // Taylor polynomial has a visible truncation error at each large step.
    for coordinate in [LocalCoordinate::Identity, LocalCoordinate::BalancedMobius] {
        let p = Precision::decimal(70)?;
        let mut opts = options(1, 8, 70);
        opts.local_coordinate = coordinate;
        let source = rational();
        let ordinary = DifferentialSystem {
            variable: source.variable,
            matrix: source.matrices[0].clone(),
        }
        .compile(p, &Default::default())?
        .transport(
            &BoundaryData {
                point: p.zero(),
                values: vec![p.i(1)],
            },
            &[p.i(1)],
            &opts,
            &RunContext::default(),
        )?;
        assert!(ordinary.diagnostics.pade_steps > 0);
        assert!(p.close(
            &ordinary.values[0],
            &p.rational(&Rational::from((1, 2))),
            45
        ));
        let result = rational().compile(p, &Default::default())?.transport(
            &EpsilonBoundary {
                point: p.zero(),
                leading: 0,
                coefficients: vec![vec![p.i(1)]],
            },
            &[p.i(1)],
            &opts,
            &RunContext::default(),
            true,
        )?;
        assert!(
            result.diagnostics.pade_steps > 0,
            "{:?}",
            result.diagnostics
        );
        assert!(p.close(
            &result.coefficients[0][0],
            &p.rational(&Rational::from((1, 2))),
            45
        ));
        for (index, segment) in result.segments.iter().enumerate() {
            assert!(segment.pade.is_some(), "{:?}", result.diagnostics);
            for point in [
                p.scale(&p.add(&segment.center, &segment.end), 1, 2),
                segment.end.clone(),
            ] {
                let expected = p.div(&p.i(1), &p.add(&p.i(1), &point));
                assert!(p.close(
                    &result.evaluate_segment(index, &point)?[0][0],
                    &expected,
                    45
                ));
            }
        }
    }
    Ok(())
}
#[test]
fn coupled_complex_epsilon_channels_agree_with_independent_precision_and_order() -> Result<()> {
    let x = symbol!("pade_coupled::x");
    let a = -Atom::one() / (Atom::one() + Atom::var(x));
    // Y'=(-I/(1+x)+epsilon*N)Y, N²=0 and N01=i.
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![
            vec![vec![a.clone(), Atom::zero()], vec![Atom::zero(), a]],
            vec![
                vec![Atom::zero(), Atom::i()],
                vec![Atom::zero(), Atom::zero()],
            ],
        ],
    };
    let target = Complex::new(Rational::from((1, 5)), Rational::from((1, 7)));
    let exact_target = Atom::num(target);
    let result = transport_epsilon(
        &system,
        |p| {
            Ok(EpsilonBoundary {
                point: p.zero(),
                leading: -1,
                coefficients: vec![vec![p.i(2), p.i(3)], vec![p.zero(), p.zero()]],
            })
        },
        &[exact_target],
        &options(2, 12, 70),
        &RunContext::default(),
        true,
    )?;
    let p = Precision {
        bits: result.diagnostics.working_bits,
    };
    let f = p.div(&p.i(1), &p.add(&p.i(1), &result.point));
    let expected = [
        p.scale(&f, 2, 1),
        p.scale(&f, 3, 1),
        p.mul(&p.complex(0, 3), &p.mul(&result.point, &f)),
        p.zero(),
    ];
    for (a, b) in result.coefficients.iter().flatten().zip(expected) {
        assert!(p.close(a, &b, 30), "{a} != {b}");
    }
    assert_eq!(result.verified_digits, Some(20));
    assert!(
        result.diagnostics.pade_steps > 0,
        "{:?}",
        result.diagnostics
    );
    assert!(!result.checkpoints.is_empty());
    Ok(())
}
#[test]
fn resource_fallback_retains_taylor_and_cancellation_commits_nothing() -> Result<()> {
    let p = Precision::decimal(60)?;
    let compiled = rational().compile(p, &Default::default())?;
    let boundary = EpsilonBoundary {
        point: p.zero(),
        leading: 0,
        coefficients: vec![vec![p.i(1)]],
    };
    let mut opts = options(4, 64, 60);
    opts.pade.as_mut().unwrap().max_work = 2;
    let result = compiled.transport(&boundary, &[p.i(1)], &opts, &RunContext::default(), true)?;
    assert_eq!(result.diagnostics.pade_steps, 0);
    assert!(result.diagnostics.pade_fallbacks > 0);
    assert!(
        result
            .diagnostics
            .last_pade_fallback
            .unwrap()
            .contains("work estimate")
    );
    assert!(result.segments.iter().all(|s| s.pade.is_none()));
    assert!(p.close(
        &result.coefficients[0][0],
        &p.rational(&Rational::from((1, 2))),
        30
    ));
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        compiled.transport(&boundary, &[p.i(1)], &opts, &context, true),
        Err(Error::Cancelled)
    ));
    assert_eq!(boundary.point, p.zero());
    assert_eq!(boundary.coefficients, vec![vec![p.i(1)]]);
    Ok(())
}
#[test]
fn cancellation_counterexample_keeps_fresh_precision_refinement() -> Result<()> {
    let x = symbol!("pade_cancellation::x");
    let system = EpsilonSystem {
        variable: x,
        matrices: vec![vec![
            vec![Atom::zero(), Atom::one() - Atom::num(10).pow(100)],
            vec![Atom::zero(), Atom::zero()],
        ]],
    };
    let calls = std::cell::RefCell::new(Vec::new());
    let result = transport_epsilon(
        &system,
        |p| {
            calls.borrow_mut().push(p.bits);
            Ok(EpsilonBoundary {
                point: p.zero(),
                leading: 0,
                coefficients: vec![vec![
                    p.eval(&Atom::num(10).pow(100), &Default::default())?,
                    p.i(1),
                ]],
            })
        },
        &[Atom::one()],
        &options(1, 16, 60),
        &RunContext::default(),
        true,
    )?;
    let p = Precision {
        bits: result.diagnostics.working_bits,
    };
    assert!(p.close(&result.coefficients[0][0], &p.i(1), 20));
    assert_eq!(result.verified_digits, Some(20));
    let calls = calls.into_inner();
    assert!(calls.len() >= 3);
    assert!(calls[1] > Precision::decimal(120)?.bits);
    Ok(())
}

#[test]
fn progressive_cache_restart_samples_the_accepted_rational_segments() -> Result<()> {
    use std::collections::BTreeMap;
    use symbolica_amflow::{kinematics::KinematicSystem, transport_cache::*};
    let s = symbol!("pade_cache::s");
    let epsilon = symbol!("pade_cache::eps");
    let flow = RustFlow::new(
        KinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(
                s,
                vec![vec![-Atom::one() / (Atom::one() + Atom::var(s))]],
            )]),
        },
        &[parse!("pade_cache::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "regular rational analytic germ",
    )?;
    let point = |n| BTreeMap::from([(s, Atom::num(n))]);
    let p = Precision::decimal(80)?;
    let range = EpsilonRange::new(0, 0)?;
    let mut bank = RustFlowCache::default();
    bank.insert(CachedBoundary {
        identity: flow.identity().clone(),
        point: CachedPoint::Exact(point(0)),
        kind: PointKind::Physical,
        range,
        coefficients: vec![vec![p.i(1)]],
        accuracy: BoundaryAccuracy::supplied(
            40,
            p.bits,
            vec![vec![p.tolerance(60)]],
            "exact analytic unit boundary",
        )?,
    })?;
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let opts = options(1, 8, 70);
    let result = flow.evaluate_to(
        &mut bank,
        &point(1),
        range,
        &opts,
        &RunContext::default(),
        &policy,
    )?;
    assert!(result.transport.as_ref().unwrap().diagnostics.pade_steps > 0);
    assert!(result.inserted_points > 1);
    for boundary in bank.entries() {
        let coordinates = boundary.point.restart_coordinates()?;
        let at = p.eval(&coordinates[&s], &Default::default())?;
        assert!(p.close(
            &boundary.coefficients[0][0],
            &p.div(&p.i(1), &p.add(&p.i(1), &at)),
            20
        ));
    }
    let dir = std::env::temp_dir().join(format!("pade-restart-{}", std::process::id()));
    bank.save(&dir)?;
    let mut restored = RustFlowCache::load(&dir)?;
    std::fs::remove_dir_all(dir)?;
    let hit = flow.evaluate_to(
        &mut restored,
        &point(1),
        range,
        &opts,
        &RunContext::default(),
        &policy,
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.inserted_points, 0);
    let before = restored.len();
    let context = RunContext::default();
    context.cancellation.cancel();
    assert!(matches!(
        flow.evaluate_to(&mut restored, &point(2), range, &opts, &context, &policy),
        Err(Error::Cancelled)
    ));
    assert_eq!(restored.len(), before);
    let next = flow.evaluate_to(
        &mut restored,
        &point(2),
        range,
        &opts,
        &RunContext::default(),
        &policy,
    )?;
    assert_eq!(next.starting_point.restart_coordinates()?, point(1));
    assert!(p.close(
        &next.boundary.coefficients[0][0],
        &p.rational(&Rational::from((1, 3))),
        20
    ));
    Ok(())
}

#[test]
fn registered_roots_explicitly_fall_back_without_changing_sheet() -> Result<()> {
    use std::collections::BTreeMap;
    use symbolica_amflow::algebraic::{AlgebraicSystem, RootSeed, SquareRoot};
    let x = symbol!("pade_roots::x");
    let r = symbol!("pade_roots::r");
    let system = AlgebraicSystem::ordinary(
        DifferentialSystem {
            variable: x,
            matrix: vec![vec![Atom::one() / Atom::var(r)]],
        },
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(x),
        }],
    );
    let p = Precision::decimal(60)?;
    let compiled = system.compile(p)?;
    let result = compiled.transport_ordinary(
        &BoundaryData {
            point: p.i(1),
            values: vec![p.i(1)],
        },
        &[
            p.complex(1, 1),
            p.complex(-1, 1),
            p.complex(-1, -1),
            p.complex(1, -1),
            p.i(1),
        ],
        &BTreeMap::from([(r, RootSeed::Principal)]),
        &options(4, 64, 60),
        &RunContext::default(),
        true,
    )?;
    assert_eq!(result.solution.diagnostics.pade_steps, 0);
    assert!(result.solution.diagnostics.pade_fallbacks > 0);
    assert!(
        result
            .solution
            .diagnostics
            .last_pade_fallback
            .unwrap()
            .contains("Taylor candidates only")
    );
    assert!(result.solution.segments.iter().all(|s| s.pade.is_none()));
    assert!(p.close(&result.branches.roots[0].1, &p.i(-1), 30));
    assert!(p.close(&result.solution.coefficients[0][0], &p.exp(&p.i(-4)), 30));
    Ok(())
}
