//! Capacity probe: padding is applied after physical source generation.
//! This verifies an embedding using public guarded APIs; production dispatch
//! and the physical weighted measure remain unchanged.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::guarded::{
    GuardedContext, GuardedDiscoveryOptions, GuardedIdentity, GuardedIdentityTerm,
    GuardedMeasureIdentity, GuardedReductionProgram, IndexBounds, IndexDomain, IndexRole,
};
use symbolica_amflow::finite_density::measure::WeightedMeasure;

const PHYSICAL: usize = 2;
const CAPACITY: usize = 4;

fn identity(capacity: usize) -> GuardedMeasureIdentity {
    GuardedMeasureIdentity {
        measure: format!(
            "integral dx x^(-a) H_b(1-x); physical_arity=2; storage_capacity={capacity}; embedding=trailing-zero-v1"
        ),
        support: "b=1 is delta(1-x); physical moments have compact point support; dummy axes exactly zero".into(),
        orientation: "h=1-x; standard real delta orientation".into(),
        normalization: "ordinary real dx; no AMF or oracle value".into(),
        branch: "polynomial test functions on the surface; theta bulk is not zero".into(),
        deformation: "none; algebraic capacity probe only".into(),
    }
}

fn derived_sources() -> Vec<GuardedIdentity<PHYSICAL>> {
    let x = parse!("capacity_probe_x");
    let measure = WeightedMeasure::new(
        vec![x.clone()],
        [x.clone(), Atom::one() - x],
        [IndexRole::Ordinary, IndexRole::Occupation],
    )
    .unwrap();
    // Derive (1-x)*delta(1-x)=0 using the actual physical owner. No source
    // coefficient or replacement source equation is authored by the probe.
    let sources: Vec<_> = measure
        .multiplication_sources()
        .unwrap()
        .into_iter()
        .filter(|source| source.domain.bounds()[1] == IndexBounds::fixed(1))
        .collect();
    assert_eq!(sources.len(), 1);
    sources
}

fn pad_domain(domain: &IndexDomain<PHYSICAL>) -> IndexDomain<CAPACITY> {
    IndexDomain::new(std::array::from_fn(|axis| {
        if axis < PHYSICAL {
            domain.bounds()[axis]
        } else {
            IndexBounds::fixed(0)
        }
    }))
    .unwrap()
}

fn pad_source(source: GuardedIdentity<PHYSICAL>) -> GuardedIdentity<CAPACITY> {
    GuardedIdentity {
        id: source.id,
        terms: source
            .terms
            .into_iter()
            .map(|term| GuardedIdentityTerm {
                shift: std::array::from_fn(|axis| term.shift.get(axis).copied().unwrap_or(0)),
                coefficient: term.coefficient,
            })
            .collect(),
        domain: pad_domain(&source.domain),
        nonzero_conditions: source.nonzero_conditions,
    }
}

fn project(point: &[i64; CAPACITY]) -> Result<[i64; PHYSICAL], &'static str> {
    if point[PHYSICAL..].iter().any(|&value| value != 0) {
        return Err("nonzero storage tail is not a physical integral");
    }
    Ok(std::array::from_fn(|axis| point[axis]))
}

fn checked_reduce(
    program: &GuardedReductionProgram<CAPACITY>,
    target: [i64; CAPACITY],
    dummy_symbols: &[Symbol],
) -> Result<BTreeMap<[i64; PHYSICAL], Atom>, &'static str> {
    project(&target)?;
    let reduced = program
        .reduce(target, Default::default())
        .map_err(|_| "native reduction error")?;
    for term in &reduced.unresolved {
        project(&term.integral)?;
    }
    if !reduced.unresolved.is_empty() {
        return Err("unresolved physical target");
    }
    let mut out = BTreeMap::new();
    for (point, coefficient) in reduced.terms {
        let point = project(&point)?;
        if dummy_symbols
            .iter()
            .any(|&symbol| !coefficient.derivative(symbol).is_zero())
        {
            return Err("coefficient depends on a dummy index");
        }
        assert!(out.insert(point, coefficient).is_none());
    }
    for condition in reduced.nonzero_conditions {
        if dummy_symbols
            .iter()
            .any(|&symbol| !condition.derivative(symbol).is_zero())
        {
            return Err("condition depends on a dummy index");
        }
    }
    Ok(out)
}

fn padded_context(
    indices: [Symbol; CAPACITY],
    measure: GuardedMeasureIdentity,
) -> GuardedContext<CAPACITY> {
    let sources = derived_sources()
        .into_iter()
        .map(pad_source)
        .collect::<Vec<_>>();
    for source in &sources {
        for axis in PHYSICAL..CAPACITY {
            assert_eq!(source.domain.bounds()[axis], IndexBounds::fixed(0));
            assert!(source.terms.iter().all(|term| term.shift[axis] == 0));
            assert!(
                source
                    .terms
                    .iter()
                    .all(|term| term.coefficient.derivative(indices[axis]).is_zero())
            );
        }
    }
    GuardedContext::new(
        measure,
        [
            IndexRole::Ordinary,
            IndexRole::Occupation,
            IndexRole::Ordinary,
            IndexRole::Ordinary,
        ],
        indices,
        vec![],
        sources,
    )
    .unwrap()
}

#[test]
fn derived_unpadded_surface_sources_reduce_and_replay_in_zero_tail_storage() {
    let indices = [
        symbol!("capacity_probe::a"),
        symbol!("capacity_probe::b"),
        symbol!("capacity_probe::dummy_2"),
        symbol!("capacity_probe::dummy_3"),
    ];
    let original = GuardedContext::new(
        identity(PHYSICAL),
        [IndexRole::Ordinary, IndexRole::Occupation],
        [indices[0], indices[1]],
        vec![],
        derived_sources(),
    )
    .unwrap();
    let padded = padded_context(indices, identity(CAPACITY));
    let domain = IndexDomain::new([
        IndexBounds::new(None, Some(-1)).unwrap(),
        IndexBounds::fixed(1),
    ])
    .unwrap();
    let options = GuardedDiscoveryOptions {
        max_depth: 2,
        max_domains: 32,
        sample_seed: 0,
    };
    let direct = original
        .discover(vec![domain.clone()], [[0, 1]], options)
        .unwrap();
    let stored = padded
        .discover(vec![pad_domain(&domain)], [[0, 1, 0, 0]], options)
        .unwrap();
    assert!(direct.unresolved.is_empty(), "{:?}", direct.unresolved);
    assert!(stored.unresolved.is_empty(), "{:?}", stored.unresolved);
    for rule in stored.program.native().rules() {
        // Discovery already certifies replay. Explicit replay also checks that
        // the program remains bound to the unmodified padded source corpus.
        padded.sources().replay_rule(rule).unwrap();
        for axis in PHYSICAL..CAPACITY {
            assert_eq!(rule.domain().bounds()[axis], IndexBounds::fixed(0));
            assert!(!rule.candidate().target[axis].is_symbolic());
            assert_eq!(rule.candidate().target[axis].value(), 0);
            for term in &rule.candidate().rhs {
                assert!(!term.integral[axis].is_symbolic());
                assert_eq!(term.integral[axis].value(), 0);
            }
        }
    }
    let bytes = stored.program.encode(Default::default()).unwrap();
    let decoded = padded.decode(&bytes, Default::default()).unwrap();
    for moment in [1, 2, 5] {
        let target = [-moment, 1];
        let expected = direct.program.reduce(target, Default::default()).unwrap();
        assert!(expected.unresolved.is_empty());
        assert_eq!(expected.terms, BTreeMap::from([([0, 1], Atom::one())]));
        let padded_target = [-moment, 1, 0, 0];
        let projected = checked_reduce(&stored.program, padded_target, &indices[2..]).unwrap();
        let restored = checked_reduce(&decoded, padded_target, &indices[2..]).unwrap();
        assert_eq!(projected, expected.terms);
        assert_eq!(restored, expected.terms);
    }
    // A surface multiplication rule cannot be applied to the bulk theta.
    assert!(checked_reduce(&decoded, [-1, 0, 0, 0], &indices[2..]).is_err());
    // Reject invalid frontend labels rather than calling them zero integrals.
    assert!(checked_reduce(&decoded, [-1, 1, 1, 0], &indices[2..]).is_err());
    assert!(project(&[0, 1, 0, -1]).is_err());
    // The same projection validates candidate terminals before discovery.
    assert!(project(&[0, 1, 1, 0]).is_err());
    let mut changed = identity(CAPACITY);
    changed
        .measure
        .push_str("; incompatible embedding version=2");
    let foreign = padded_context(indices, changed);
    assert!(foreign.decode(&bytes, Default::default()).is_err());
}
