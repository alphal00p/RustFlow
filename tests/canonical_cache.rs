use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::algebraic::{AlgebraicKinematicSystem, CanonicalAlgebraicSystem, SquareRoot};
use symbolica_amflow::transport_cache::*;
use symbolica_amflow::{Error, FlowOptions, Precision, Prescription, Result, RunContext, RustFlow};

fn names() -> (Symbol, Symbol, Symbol) {
    (
        symbol!("canonical_cache::eps"),
        symbol!("canonical_cache::s"),
        symbol!("canonical_cache::r"),
    )
}
fn canonical() -> Result<CanonicalAlgebraicSystem> {
    let (eps, s, r) = names();
    CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::var(s) + Atom::var(r) + 1],
        &[vec![vec![Atom::one()]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    )
}
fn flow() -> Result<RustFlow<CanonicalAlgebraicSystem>> {
    RustFlow::new_canonical(
        canonical()?,
        &[Atom::var(symbol!("canonical_cache::I"))],
        &Atom::one(),
        Prescription::PlusI0,
        "positive real s; no additional monodromy",
    )
}
fn germ(sheet: RootSheet) -> RootGerm {
    RootGerm {
        sheets: BTreeMap::from([(names().2, sheet)]),
    }
}
fn point(value: i64) -> BTreeMap<Symbol, Atom> {
    BTreeMap::from([(names().1, Atom::num(value))])
}
fn seed(identity: &BoundaryIdentity, value: i64, digits: u32) -> Result<CachedBoundary> {
    let p = Precision::decimal(80)?;
    let point = CachedPoint::Exact(point(value));
    let point = if identity.roots().is_empty() {
        point
    } else {
        point.with_root_germ(germ(RootSheet::Principal))?
    };
    Ok(CachedBoundary {
        identity: identity.clone(),
        point,
        kind: PointKind::Physical,
        range: EpsilonRange::new(0, 2)?,
        coefficients: vec![vec![p.i(1)], vec![p.zero()], vec![p.zero()]],
        accuracy: BoundaryAccuracy::supplied(
            digits,
            p.bits,
            vec![vec![p.real(0)]; 3],
            "exact independent normalization Y(source)=1",
        )?,
    })
}
fn options() -> FlowOptions {
    FlowOptions {
        digits: 20,
        guard_digits: 30,
        series_order: 64,
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
fn canonical_cache_reuses_restarts_and_agrees_with_refined_dense_transport() -> Result<()> {
    let flow = flow()?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(flow.identity(), 1, 35)?)?;
    let range = EpsilonRange::new(0, 2)?;
    let p = Precision::decimal(80)?;
    let germ = germ(RootSheet::Principal);
    let first = flow.evaluate_to(
        &mut bank,
        &point(4),
        Some(&germ),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(first.transport.is_some());
    let result = flow.evaluate_to(
        &mut bank,
        &point(9),
        Some(&germ),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert_eq!(result.starting_point.restart_coordinates()?, point(4));
    let logarithm = p.log(&p.div(&p.i(13), &p.i(3)));
    for (actual, expected) in result.boundary.coefficients.iter().zip([
        p.i(1),
        logarithm.clone(),
        p.div(&p.mul(&logarithm, &logarithm), &p.i(2)),
    ]) {
        assert!(p.close(&actual[0], &expected, 20));
    }
    let canonical = canonical()?;
    let dense = AlgebraicKinematicSystem::canonical_dlog(
        canonical.epsilon(),
        canonical.variables(),
        canonical.letters(),
        canonical.constant_matrices(),
        canonical.roots().to_vec(),
    )?;
    let dense = RustFlow::new_algebraic(
        dense,
        &[Atom::var(symbol!("canonical_cache::I"))],
        &Atom::one(),
        Prescription::PlusI0,
        "positive real s; no additional monodromy",
    )?;
    assert_ne!(dense.identity().key(), flow.identity().key());
    let mut independent = RustFlowCache::default();
    independent.insert(seed(dense.identity(), 1, 35)?)?;
    let refined = FlowOptions {
        guard_digits: 50,
        series_order: 96,
        ..options()
    };
    let dense_result = dense.evaluate_to(
        &mut independent,
        &point(9),
        &germ,
        range,
        &refined,
        &RunContext::default(),
        &policy(),
    )?;
    for (a, b) in result
        .boundary
        .coefficients
        .iter()
        .flatten()
        .zip(dense_result.boundary.coefficients.iter().flatten())
    {
        assert!(p.close(a, b, 20));
    }
    // One physical bank can retain both explicitly different representations.
    bank.insert(dense_result.boundary)?;
    let directory = std::env::temp_dir().join(format!("canonical-cache-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    bank.save(&directory)?;
    let mut restored = RustFlowCache::load(&directory)?;
    assert_eq!(bank.len(), restored.len());
    let exact = flow.evaluate_to(
        &mut restored,
        &point(9),
        Some(&germ),
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    assert!(exact.transport.is_none());
    assert_eq!(exact.boundary.identity.key(), flow.identity().key());
    assert!(matches!(
        flow.evaluate_to(
            &mut restored,
            &point(9),
            Some(&self::germ(RootSheet::Opposite)),
            range,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::IncompleteReduction(_))
    ));
    assert!(matches!(
        flow.evaluate_to(
            &mut restored,
            &point(9),
            None,
            range,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::InvalidInput(_))
    ));
    std::fs::remove_dir_all(directory)?;
    Ok(())
}

#[test]
fn canonical_identity_preserves_order_and_never_collides_with_dense_zero_matrices() -> Result<()> {
    let (eps, s, r) = names();
    let roots = vec![SquareRoot {
        symbol: r,
        radicand: Atom::var(s),
    }];
    let letters = [Atom::var(s) + Atom::var(r) + 1, Atom::var(s)];
    let matrices = [vec![vec![Atom::one()]], vec![vec![Atom::num(2)]]];
    let a = CanonicalAlgebraicSystem::new(eps, &[s], &letters, &matrices, roots.clone())?;
    let b = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[letters[1].clone(), letters[0].clone()],
        &[matrices[1].clone(), matrices[0].clone()],
        roots.clone(),
    )?;
    let build = |system| {
        RustFlow::new_canonical(
            system,
            &[Atom::one()],
            &Atom::one(),
            Prescription::PlusI0,
            "same exact sheet",
        )
    };
    assert_ne!(build(a)?.identity().key(), build(b)?.identity().key());
    let zero = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::var(s)],
        &[vec![vec![Atom::new()]]],
        roots.clone(),
    )?;
    let dense = RustFlow::new_algebraic(
        AlgebraicKinematicSystem {
            epsilon: eps,
            derivatives: BTreeMap::from([(s, vec![vec![Atom::new()]])]),
            roots,
        },
        &[Atom::one()],
        &Atom::one(),
        Prescription::PlusI0,
        "same exact sheet",
    )?;
    assert_ne!(build(zero)?.identity().key(), dense.identity().key());
    Ok(())
}

#[test]
fn rootless_canonical_transport_uses_the_same_cache_without_a_germ() -> Result<()> {
    let (eps, s, _) = names();
    let canonical = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[Atom::var(s)],
        &[vec![vec![Atom::one()]]],
        vec![],
    )?;
    let flow = RustFlow::new_canonical(
        canonical,
        &[Atom::one()],
        &Atom::one(),
        Prescription::PlusI0,
        "positive real",
    )?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(flow.identity(), 1, 35)?)?;
    let range = EpsilonRange::new(0, 2)?;
    let result = flow.evaluate_to(
        &mut bank,
        &point(4),
        None,
        range,
        &options(),
        &RunContext::default(),
        &policy(),
    )?;
    let p = Precision::decimal(80)?;
    let logarithm = p.log(&p.i(4));
    assert!(result.boundary.point.root_germ().is_none());
    assert!(p.close(&result.boundary.coefficients[1][0], &logarithm, 20));
    assert!(p.close(
        &result.boundary.coefficients[2][0],
        &p.div(&p.mul(&logarithm, &logarithm), &p.i(2)),
        20
    ));
    assert!(matches!(
        flow.evaluate_to(
            &mut bank,
            &point(4),
            Some(&germ(RootSheet::Principal)),
            range,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::InvalidInput(_))
    ));
    Ok(())
}

#[test]
fn canonical_canceled_letters_keep_path_holes_and_supplied_uncertainty() -> Result<()> {
    let (eps, s, r) = names();
    let letter = Atom::var(s) - 4;
    let system = CanonicalAlgebraicSystem::new(
        eps,
        &[s],
        &[letter.clone(), letter],
        &[vec![vec![Atom::one()]], vec![vec![Atom::num(-1)]]],
        vec![SquareRoot {
            symbol: r,
            radicand: Atom::var(s),
        }],
    )?;
    let canceled = RustFlow::new_canonical(
        system,
        &[Atom::one()],
        &Atom::one(),
        Prescription::PlusI0,
        "positive real with source holes",
    )?;
    let mut bank = RustFlowCache::default();
    bank.insert(seed(canceled.identity(), 1, 35)?)?;
    let range = EpsilonRange::new(0, 2)?;
    let germ = germ(RootSheet::Principal);
    assert!(matches!(
        canceled.evaluate_to(
            &mut bank,
            &point(9),
            Some(&germ),
            range,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::IncompleteReduction(_))
    ));
    assert_eq!(bank.len(), 1);
    assert!(matches!(
        canceled.evaluate_to(
            &mut bank,
            &point(4),
            Some(&germ),
            range,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::InvalidInput(_))
    ));
    let flow = flow()?;
    let mut weak = RustFlowCache::default();
    weak.insert(seed(flow.identity(), 1, 20)?)?;
    assert!(matches!(
        flow.evaluate_to(
            &mut weak,
            &point(4),
            Some(&germ),
            range,
            &options(),
            &RunContext::default(),
            &policy()
        ),
        Err(Error::Accuracy(_))
    ));
    assert_eq!(weak.len(), 1);
    Ok(())
}
