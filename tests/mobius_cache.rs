use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, SquareRoot};
use symbolica_amflow::local_coordinates::TaylorCoordinate;
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::*;

fn symbols() -> (Symbol, Symbol, Symbol) {
    (
        symbol!("mobius_bank::s"),
        symbol!("mobius_bank::eps"),
        symbol!("mobius_bank::r"),
    )
}
fn germ(sheet: RootSheet) -> RootGerm {
    RootGerm {
        sheets: BTreeMap::from([(symbols().2, sheet)]),
    }
}
fn flow(radicand: Atom, conditions: &[Atom]) -> Result<RustFlow<AlgebraicKinematicSystem>> {
    let (s, epsilon, r) = symbols();
    RustFlow::with_algebraic_conditions(
        AlgebraicKinematicSystem {
            epsilon,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::var(epsilon) / Atom::var(r)]])]),
            roots: vec![SquareRoot {
                symbol: r,
                radicand,
            }],
        },
        &[parse!("mobius_bank::I")],
        &Atom::one(),
        Prescription::PlusI0,
        "specified local root sheet in a simply connected regular real component",
        conditions,
    )
}
fn point(value: Atom, sheet: RootSheet) -> Result<CachedPoint> {
    CachedPoint::Exact(BTreeMap::from([(symbols().0, value)])).with_root_germ(germ(sheet))
}
fn boundary(
    flow: &RustFlow<AlgebraicKinematicSystem>,
    value: Atom,
    sheet: RootSheet,
    coefficients: &[i64],
    digits: u32,
) -> Result<CachedBoundary> {
    let p = Precision::decimal(80)?;
    Ok(CachedBoundary {
        identity: flow.identity().clone(),
        point: point(value, sheet)?,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, coefficients.len() as i32 - 1)?,
        coefficients: coefficients.iter().map(|&a| vec![p.i(a)]).collect(),
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.tolerance(60)]; coefficients.len()],
            "independent analytic boundary",
        )?,
    })
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 30,
        series_order: 64,
        local_coordinate: LocalCoordinate::BalancedMobius,
        ..Default::default()
    }
}
fn policy() -> impl TransportCost {
    ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    }
}

#[test]
fn mapped_physical_checkpoints_preserve_germs_and_restart_binary() -> Result<()> {
    let s = symbols().0;
    let flow = flow(Atom::var(s), &[])?;
    let mut cache = RustFlowCache::default();
    cache.insert(boundary(
        &flow,
        Atom::one(),
        RootSheet::Principal,
        &[1, 0, 0],
        35,
    )?)?;
    let range = EpsilonRange::new(0, 2)?;
    let first = flow.evaluate_to(
        &mut cache,
        &BTreeMap::from([(s, Atom::num(4))]),
        &germ(RootSheet::Principal),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(first.inserted_points > 1);
    assert!(
        first
            .transport
            .as_ref()
            .unwrap()
            .segments
            .iter()
            .any(|segment| matches!(segment.coordinate, TaylorCoordinate::Mobius(_)))
    );
    assert!(cache.entries().iter().any(|b|matches!(&b.point,CachedPoint::Algebraic {point,..} if matches!(point.as_ref(),CachedPoint::Derived{..}))));
    assert!(
        cache
            .entries()
            .iter()
            .all(|b| b.point.root_germ() == Some(&germ(RootSheet::Principal)))
    );
    let path = std::env::temp_dir().join(format!("rustflow-mobius-cache-{}", std::process::id()));
    cache.save(&path)?;
    let mut restored = RustFlowCache::load(&path)?;
    std::fs::remove_dir_all(&path)?;
    assert_eq!(cache.len(), restored.len());
    let second = flow.evaluate_to(
        &mut restored,
        &BTreeMap::from([(s, Atom::num(9))]),
        &germ(RootSheet::Principal),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!((second.starting_point.restart_coordinates()?[&s].clone() - 4).is_zero());
    let p = Precision::decimal(80)?;
    // exp(2 epsilon (sqrt(s)-1)) gives [1,4,8] at s=9.
    for (actual, expected) in second.boundary.coefficients.iter().flatten().zip([1, 4, 8]) {
        assert!(p.close(actual, &p.i(expected), 30));
    }
    let hit = flow.evaluate_to(
        &mut restored,
        &BTreeMap::from([(s, Atom::num(9))]),
        &germ(RootSheet::Principal),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(hit.transport.is_none());
    assert_eq!(hit.inserted_points, 0);
    Ok(())
}
