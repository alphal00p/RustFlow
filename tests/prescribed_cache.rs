use std::{collections::BTreeMap, sync::Arc};
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, CanonicalAlgebraicSystem, SquareRoot};
use symbolica_amflow::contour::PolynomialPrescription;
use symbolica_amflow::kinematics::KinematicSystem;
use symbolica_amflow::physical_transport::PhysicalRoute;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{Error, FlowOptions, Precision, Prescription, Result, RunContext, RustFlow};

fn names() -> (Symbol, Symbol, Symbol) {
    (
        symbol!("prescribed_cache::s"),
        symbol!("prescribed_cache::eps"),
        symbol!("prescribed_cache::r"),
    )
}
fn point(n: i64) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(names().0, Atom::num(n))])
}
fn germ(sheet: RootSheet) -> RootGerm {
    RootGerm {
        sheets: BTreeMap::from([(names().2, sheet)]),
    }
}
fn convention(side: Prescription, domain: &str) -> PhysicalContinuation {
    PhysicalContinuation {
        prescriptions: vec![PolynomialPrescription {
            polynomial: Atom::var(names().0),
            prescription: side,
        }],
        unprescribed_side: side,
        domain: domain.into(),
    }
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 30,
        series_order: 48,
        ..Default::default()
    }
}
fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}
fn admit(_: &CachedBoundary, _: &CachedPoint, _: &PhysicalRoute) -> Result<bool> {
    Ok(true)
}
fn seed<S>(flow: &RustFlow<S>, at: i64, sheet: Option<RootSheet>) -> Result<CachedBoundary> {
    let p = Precision::decimal(100)?;
    let point = CachedPoint::Exact(point(at));
    let point = if let Some(sheet) = sheet {
        point.with_root_germ(germ(sheet))?
    } else {
        point
    };
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 2)?,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            70,
            p.bits,
            vec![vec![p.real(0)]; 3],
            "exact normalization Y(source)=1",
        )?,
    })
}
fn root_flow() -> Result<RustFlow<AlgebraicKinematicSystem>> {
    let (s, e, r) = names();
    RustFlow::new_algebraic(
        AlgebraicKinematicSystem {
            epsilon: e,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::var(e) / Atom::var(r)]])]),
            roots: vec![SquareRoot {
                symbol: r,
                radicand: Atom::var(s),
            }],
        },
        &[parse!("prescribed_cache::Y")],
        &Atom::one(),
        Prescription::MinusI0,
        "sqrt branch normalized at positive s",
    )?
    .with_prescribed_continuation(convention(
        Prescription::MinusI0,
        "lower physical s half-plane; no extra winding",
    ))
}
fn cache_bytes(bank: &RustFlowCache, label: &str) -> Result<Vec<u8>> {
    let dir = std::env::temp_dir().join(format!("prescribed-bank-{}-{label}", std::process::id()));
    bank.save(&dir)?;
    let file = std::fs::read_dir(&dir)?.next().unwrap()?.path();
    let bytes = std::fs::read(file)?;
    std::fs::remove_dir_all(dir)?;
    Ok(bytes)
}

#[test]
fn root_threshold_progressive_points_have_actual_germs_and_restart() -> Result<()> {
    let flow = root_flow()?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(&flow, 1, Some(RootSheet::Principal))?)?;
    let result = flow.evaluate_prescribed_to(
        &mut bank,
        &point(-1),
        &germ(RootSheet::Opposite),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    let p = Precision::decimal(90)?;
    for (row, expected) in
        result
            .boundary
            .coefficients
            .iter()
            .zip([p.i(1), p.complex(-2, -2), p.complex(0, 4)])
    {
        assert!(
            p.close(&row[0], &expected, 20),
            "{} vs {}",
            row[0],
            expected
        );
    }
    assert!(result.inserted_points > 2);
    let mut before = false;
    let mut after = false;
    for entry in bank.entries() {
        assert_eq!(entry.kind, PointKind::Physical);
        let s = p.eval(
            &entry.point.restart_coordinates()?[&names().0],
            &Default::default(),
        )?;
        assert_eq!(s.im, p.real(0));
        let expected = if s.re > p.real(0) {
            before = true;
            RootSheet::Principal
        } else {
            after = true;
            RootSheet::Opposite
        };
        assert_eq!(entry.point.root_germ(), Some(&germ(expected)));
    }
    assert!(before && after);
    let dir = std::env::temp_dir().join(format!("prescribed-restart-{}", std::process::id()));
    bank.save(&dir)?;
    let mut restored = RustFlowCache::load(&dir)?;
    std::fs::remove_dir_all(dir)?;
    let second = flow.evaluate_prescribed_to(
        &mut restored,
        &point(-4),
        &germ(RootSheet::Opposite),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    assert_eq!(second.starting_point.restart_coordinates()?, point(-1));
    for (row, expected) in
        second
            .boundary
            .coefficients
            .iter()
            .zip([p.i(1), p.complex(-2, -4), p.complex(-6, 8)])
    {
        assert!(p.close(&row[0], &expected, 20));
    }
    let hit = flow.evaluate_prescribed_to(
        &mut restored,
        &point(-4),
        &germ(RootSheet::Opposite),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    assert!(hit.transport.is_none());
    let mut independent = RustFlowCache::default();
    independent.insert(seed(&flow, 1, Some(RootSheet::Principal))?)?;
    let refined = FlowOptions {
        guard_digits: 50,
        series_order: 80,
        ..options()
    };
    let direct = flow.evaluate_prescribed_to(
        &mut independent,
        &point(-4),
        &germ(RootSheet::Opposite),
        EpsilonRange::new(0, 2)?,
        &refined,
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    for (a, b) in second
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(direct.boundary.coefficients.iter().flatten())
    {
        assert!(p.close(a, b, 20));
    }
    Ok(())
}

#[test]
fn log_crossing_uses_typed_identity_and_cannot_alias_a_different_monodromy_domain() -> Result<()> {
    let (s, e, _) = names();
    let base = || {
        RustFlow::new(
            KinematicSystem {
                epsilon: e,
                derivatives: BTreeMap::from([(s, vec![vec![Atom::var(e) / Atom::var(s)]])]),
            },
            &[parse!("prescribed_cache::Y")],
            &Atom::one(),
            Prescription::PlusI0,
            "principal logarithm at s<0+i0",
        )
    };
    let ordinary = base()?;
    let upper = base()?.with_prescribed_continuation(convention(
        Prescription::PlusI0,
        "no additional logarithmic winding",
    ))?;
    let lower = base()?.with_prescribed_continuation(convention(
        Prescription::MinusI0,
        "no additional logarithmic winding",
    ))?;
    let wound = base()?.with_prescribed_continuation(convention(
        Prescription::PlusI0,
        "one additional logarithmic winding",
    ))?;
    assert_ne!(upper.identity().key(), ordinary.identity().key());
    assert_ne!(upper.identity().key(), lower.identity().key());
    assert_ne!(upper.identity().key(), wound.identity().key());
    let mut bank = RustFlowCache::default();
    bank.insert(seed(&upper, -1, None)?)?;
    let result = upper.evaluate_prescribed_to(
        &mut bank,
        &point(1),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    let p = Precision::decimal(90)?;
    let log = p.complex(0, -1);
    let log = p.mul(
        &log,
        &symbolica_amflow::ComplexFloat::new(p.real(1).pi(), p.real(0)),
    );
    for (row, expected) in result.boundary.coefficients.iter().zip([
        p.i(1),
        log.clone(),
        p.scale(&p.mul(&log, &log), 1, 2),
    ]) {
        assert!(p.close(&row[0], &expected, 20));
    }
    assert!(matches!(
        upper.evaluate_to(
            &mut bank,
            &point(1),
            EpsilonRange::new(0, 2)?,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::InvalidInput(_))
    ));
    assert!(matches!(
        wound.evaluate_prescribed_to(
            &mut bank,
            &point(1),
            EpsilonRange::new(0, 2)?,
            &options(),
            &RunContext::default(),
            &policy(),
            &admit
        ),
        Err(Error::IncompleteReduction(_))
    ));
    assert!(matches!(
        lower.evaluate_prescribed_to(
            &mut bank,
            &point(1),
            EpsilonRange::new(0, 2)?,
            &options(),
            &RunContext::default(),
            &policy(),
            &admit
        ),
        Err(Error::IncompleteReduction(_))
    ));
    Ok(())
}

#[test]
fn wrong_sheet_denied_homotopy_and_cancellation_leave_the_bank_unchanged() -> Result<()> {
    let flow = root_flow()?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(&flow, 1, Some(RootSheet::Principal))?)?;
    let old = cache_bytes(&bank, "initial")?;
    let wrong = flow.evaluate_prescribed_to(
        &mut bank,
        &point(-1),
        &germ(RootSheet::Principal),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    );
    assert!(
        matches!(wrong,Err(Error::InvalidInput(message))if message.contains("different root germ"))
    );
    assert_eq!(cache_bytes(&bank, "wrong")?, old);
    let deny = |_: &CachedBoundary, _: &CachedPoint, _: &PhysicalRoute| Ok(false);
    assert!(matches!(
        flow.evaluate_prescribed_to(
            &mut bank,
            &point(-1),
            &germ(RootSheet::Opposite),
            EpsilonRange::new(0, 2)?,
            &options(),
            &RunContext::default(),
            &policy(),
            &deny
        ),
        Err(Error::IncompleteReduction(_))
    ));
    assert_eq!(cache_bytes(&bank, "deny")?, old);
    let token = symbolica_amflow::CancellationToken::default();
    let cancel = token.clone();
    let context = RunContext {
        cancellation: token,
        progress: Some(Arc::new(move |event| {
            if matches!(event, symbolica_amflow::Progress::Step { index: 2 }) {
                cancel.cancel();
            }
        })),
    };
    assert!(matches!(
        flow.evaluate_prescribed_to(
            &mut bank,
            &point(-1),
            &germ(RootSheet::Opposite),
            EpsilonRange::new(0, 2)?,
            &options(),
            &context,
            &policy(),
            &admit
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(cache_bytes(&bank, "cancel")?, old);
    Ok(())
}

#[test]
fn original_guard_crossings_and_malformed_physical_declarations_are_not_erased() -> Result<()> {
    let (s, e, r) = names();
    let flow = RustFlow::with_conditions(
        KinematicSystem {
            epsilon: e,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]])]),
        },
        &[parse!("prescribed_cache::Y")],
        &Atom::one(),
        Prescription::PlusI0,
        "guarded identity",
        &[Atom::var(s)],
    )?;
    for expression in [
        Atom::var(e),
        Atom::var(r),
        parse!("prescribed_cache::foreign"),
        Atom::var(s).pow(Atom::num((1, 2))),
    ] {
        let mut bad = convention(Prescription::PlusI0, "invalid");
        bad.prescriptions[0].polynomial = expression;
        assert!(flow.identity().with_prescribed_continuation(bad).is_err());
    }
    let flow = flow.with_prescribed_continuation(convention(
        Prescription::PlusI0,
        "preserved original guard",
    ))?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(&flow, -1, None)?)?;
    let checked = std::cell::Cell::new(false);
    let admission = |_: &CachedBoundary, _: &CachedPoint, route: &PhysicalRoute| {
        assert_eq!(route.crossings.len(), 1);
        assert!(route.waypoints.len() > 1);
        checked.set(true);
        Ok(true)
    };
    let result = flow.evaluate_prescribed_to(
        &mut bank,
        &point(1),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admission,
    )?;
    assert!(checked.get());
    assert_eq!(
        result.boundary.coefficients[0][0],
        Precision::decimal(90)?.i(1)
    );
    let old = cache_bytes(&bank, "guard-good")?;
    assert!(
        flow.evaluate_prescribed_to(
            &mut bank,
            &point(0),
            EpsilonRange::new(0, 2)?,
            &options(),
            &RunContext::default(),
            &policy(),
            &admit
        )
        .is_err()
    );
    assert_eq!(cache_bytes(&bank, "guard-zero")?, old);
    Ok(())
}

#[test]
fn canonical_registered_root_uses_the_same_prescribed_cache_path() -> Result<()> {
    let (s, e, r) = names();
    let system = CanonicalAlgebraicSystem::new(
        e,
        &[s],
        &[Atom::var(r)],
        &[vec![vec![Atom::num(2)]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    )?;
    let flow = RustFlow::new_canonical(
        system,
        &[parse!("prescribed_cache::Y")],
        &Atom::one(),
        Prescription::MinusI0,
        "root logarithm",
    )?
    .with_prescribed_continuation(convention(
        Prescription::MinusI0,
        "lower half-plane, no extra windings",
    ))?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(&flow, 1, Some(RootSheet::Principal))?)?;
    let result = flow.evaluate_prescribed_to(
        &mut bank,
        &point(-1),
        Some(&germ(RootSheet::Opposite)),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    let p = Precision::decimal(90)?;
    let z = p.mul(
        &p.complex(0, -1),
        &symbolica_amflow::ComplexFloat::new(p.real(1).pi(), p.real(0)),
    );
    for (row, expected) in
        result
            .boundary
            .coefficients
            .iter()
            .zip([p.i(1), z.clone(), p.scale(&p.mul(&z, &z), 1, 2)])
    {
        assert!(p.close(&row[0], &expected, 20));
    }
    Ok(())
}

#[test]
fn constant_opposite_sheet_candidate_cannot_mask_a_valid_farther_source() -> Result<()> {
    let flow = root_flow()?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(&flow, 1, Some(RootSheet::Principal))?)?;
    bank.insert(seed(&flow, -1, Some(RootSheet::Principal))?)?;
    let result = flow.evaluate_prescribed_to(
        &mut bank,
        &point(-1),
        &germ(RootSheet::Opposite),
        EpsilonRange::new(0, 2)?,
        &options(),
        &RunContext::default(),
        &policy(),
        &admit,
    )?;
    assert_eq!(result.starting_point.restart_coordinates()?, point(1));
    assert!(Precision::decimal(90)?.close(
        &result.boundary.coefficients[1][0],
        &Precision::decimal(90)?.complex(-2, -2),
        20
    ));
    Ok(())
}
