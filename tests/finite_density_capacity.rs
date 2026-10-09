//! Capacity probe: padding is applied after physical source generation.
//! The production adapter validates zero storage tails at source construction,
//! discovery, replay, application and shared-engine export boundaries.
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
    GuardedContext::new_with_physical_arity(
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
        PHYSICAL,
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
    assert_eq!(padded.physical_arity(), PHYSICAL);
    assert!(
        padded
            .sources()
            .measure_id()
            .ends_with("zero-tail-storage-v1:physical=2:capacity=4")
    );
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
    // Production rejects invalid labels before any native zero/terminal shortcut.
    assert!(decoded.reduce([-1, 1, 1, 0], Default::default()).is_err());
    assert!(
        padded
            .discover(vec![pad_domain(&domain)], [[0, 1, 1, 0]], options)
            .is_err()
    );
    let mut invalid_domain = *pad_domain(&domain).bounds();
    invalid_domain[2] = IndexBounds::new(None, Some(0)).unwrap();
    assert!(
        padded
            .discover(
                vec![IndexDomain::new(invalid_domain).unwrap()],
                [[0, 1, 0, 0]],
                options
            )
            .is_err()
    );
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

#[test]
fn padded_source_admission_rejects_every_tail_escape_and_binds_arity() {
    let indices = [
        symbol!("capacity_admission_a"),
        symbol!("capacity_admission_b"),
        symbol!("capacity_admission_tail2"),
        symbol!("capacity_admission_tail3"),
    ];
    let roles = [
        IndexRole::Ordinary,
        IndexRole::Occupation,
        IndexRole::Ordinary,
        IndexRole::Ordinary,
    ];
    let sources = derived_sources()
        .into_iter()
        .map(pad_source)
        .collect::<Vec<_>>();
    let create = |roles, sources, physical_arity| {
        GuardedContext::new_with_physical_arity(
            identity(CAPACITY),
            roles,
            indices,
            vec![],
            sources,
            physical_arity,
        )
    };
    assert!(create(roles, sources.clone(), 0).is_err());
    assert!(create(roles, sources.clone(), CAPACITY + 1).is_err());
    let mut bad_roles = roles;
    bad_roles[2] = IndexRole::Occupation;
    assert!(create(bad_roles, sources.clone(), PHYSICAL).is_err());
    let mut shifted = sources.clone();
    shifted[0].terms[0].shift[2] = 1;
    assert!(create(roles, shifted, PHYSICAL).is_err());
    let mut unbounded = sources.clone();
    let mut bounds = *unbounded[0].domain.bounds();
    bounds[2] = IndexBounds::new(None, Some(0)).unwrap();
    unbounded[0].domain = IndexDomain::new(bounds).unwrap();
    assert!(create(roles, unbounded, PHYSICAL).is_err());
    let mut symbolic = sources.clone();
    symbolic[0].terms[0].coefficient *= Atom::var(indices[2]);
    assert!(create(roles, symbolic, PHYSICAL).is_err());
    let mut guarded = sources.clone();
    guarded[0]
        .nonzero_conditions
        .push(Atom::var(indices[2]) + Atom::one());
    assert!(create(roles, guarded, PHYSICAL).is_err());
    let context = create(roles, sources.clone(), PHYSICAL).unwrap();
    let domain = IndexDomain::new([
        IndexBounds::new(None, Some(-1)).unwrap(),
        IndexBounds::fixed(1),
        IndexBounds::fixed(0),
        IndexBounds::fixed(0),
    ])
    .unwrap();
    let found = context
        .discover(vec![domain], [[0, 1, 0, 0]], Default::default())
        .unwrap();
    let bytes = found.program.encode(Default::default()).unwrap();
    // The same storage corpus interpreted as three physical factors is a
    // different owner even if the third coordinate happens to be frozen.
    let other_arity = create(roles, sources, PHYSICAL + 1).unwrap();
    assert!(other_arity.decode(&bytes, Default::default()).is_err());
}

fn radial_closed<const N: usize>()
-> symbolica_amflow::finite_density::reduction::WeightedReducedSystem<N> {
    use symbolica_amflow::finite_density::reduction::{
        AuxiliaryConvention, FixedShellDeformation, WeightedClosureOutcome, prepare_weighted_system,
    };
    let z = symbol!("capacity_radial_z");
    let t = symbol!("capacity_radial_t");
    let radius = symbol!("capacity_radial_R");
    let dimension = symbol!("capacity_radial_d");
    let indices = std::array::from_fn(|axis| {
        match Atom::parse(
            &format!("capacity_radial_index_{axis}"),
            "rustflow_capacity",
            Default::default(),
        )
        .unwrap()
        .as_view()
        {
            AtomView::Var(variable) => variable.get_symbol(),
            _ => unreachable!(),
        }
    });
    let measure = WeightedMeasure::<N>::from_physical(
        vec![Atom::var(z)],
        vec![
            Atom::var(z) + Atom::var(t),
            Atom::var(radius) - Atom::var(z),
        ],
        vec![IndexRole::Ordinary, IndexRole::Occupation],
    )
    .unwrap();
    assert_eq!(measure.physical_arity(), 2);
    let mut sources = measure
        .ibp(
            "radial-dilation",
            &indices,
            &[Atom::num(2) * Atom::var(z)],
            &Atom::var(dimension),
            2,
        )
        .unwrap();
    sources.extend(measure.multiplication_sources().unwrap());
    let mut metadata = identity(PHYSICAL);
    metadata.measure = "radial d-dimensional ball: dz z^(d/2-1) (z+t)^(-a) H_b(R-z)".into();
    metadata.support = "R,t>0; Re(d)>0; zero storage tail".into();
    metadata.deformation = "Euclidean z+t, fixed R-z".into();
    let context = GuardedContext::new_with_physical_arity(
        metadata,
        *measure.roles(),
        indices,
        vec![t, radius, dimension],
        sources,
        PHYSICAL,
    )
    .unwrap();
    let deformation = FixedShellDeformation::new(
        t,
        AuxiliaryConvention::EuclideanPlusT,
        *measure.roles(),
        std::array::from_fn(|axis| axis == 0),
    )
    .unwrap()
    .with_physical_arity(PHYSICAL)
    .unwrap();
    let target = std::array::from_fn(|axis| if axis == 0 { 1 } else { 0 });
    for &index_symbol in &indices {
        let mut invalid_variable = deformation.clone();
        invalid_variable.variable = index_symbol;
        assert!(
            prepare_weighted_system(
                &context,
                &[BTreeMap::from([(target, Atom::one())])],
                &invalid_variable,
                Default::default(),
                &Default::default()
            )
            .is_err()
        );
    }
    if N > PHYSICAL {
        let mut invalid = target;
        invalid[PHYSICAL] = 1;
        assert!(deformation.derivative(invalid).is_err());
        assert!(
            prepare_weighted_system(
                &context,
                &[BTreeMap::from([(invalid, Atom::new())])],
                &deformation,
                Default::default(),
                &Default::default()
            )
            .is_err()
        );
        assert!(
            deformation
                .clone()
                .with_admitted_domain(IndexDomain::for_roles(measure.roles()))
                .is_err()
        );
        let mismatched = FixedShellDeformation::new(
            t,
            AuxiliaryConvention::EuclideanPlusT,
            *measure.roles(),
            std::array::from_fn(|axis| axis == 0),
        )
        .unwrap();
        assert!(
            prepare_weighted_system(
                &context,
                &[BTreeMap::from([(target, Atom::one())])],
                &mismatched,
                Default::default(),
                &Default::default()
            )
            .is_err()
        );
        assert!(
            FixedShellDeformation::new(
                t,
                AuxiliaryConvention::EuclideanPlusT,
                *measure.roles(),
                std::array::from_fn(|axis| axis == PHYSICAL)
            )
            .unwrap()
            .with_physical_arity(PHYSICAL)
            .is_err()
        );
    }
    let outcome = prepare_weighted_system(
        &context,
        &[BTreeMap::from([(target, Atom::one())])],
        &deformation,
        Default::default(),
        &Default::default(),
    )
    .unwrap();
    let WeightedClosureOutcome::Closed(closed) = outcome else {
        panic!("native radial closure failed: {outcome:?}")
    };
    assert_eq!(closed.physical_arity(), PHYSICAL);
    for basis in &closed.reduced.basis {
        assert_eq!(basis.0.len(), PHYSICAL);
    }
    for target in &closed.reduced.targets {
        for index in target.keys() {
            assert_eq!(index.0.len(), PHYSICAL);
        }
    }
    for (index, terms) in &closed.reduced.candidates {
        assert_eq!(index.0.len(), PHYSICAL);
        for index in terms.keys() {
            assert_eq!(index.0.len(), PHYSICAL);
        }
    }
    let decoded = context
        .decode(
            &closed.program.encode(Default::default()).unwrap(),
            Default::default(),
        )
        .unwrap();
    for basis in &closed.reduced.basis {
        let index =
            std::array::from_fn(|axis| basis.0.get(axis).copied().map(i64::from).unwrap_or(0));
        for derivative in deformation.derivative(index).unwrap().keys() {
            assert!(
                decoded
                    .reduce(*derivative, Default::default())
                    .unwrap()
                    .unresolved
                    .is_empty()
            );
        }
    }
    closed
}

#[test]
fn physical_basis_targets_candidates_and_connection_are_identical_after_padding() {
    let direct = radial_closed::<2>();
    let padded = radial_closed::<4>();
    assert_eq!(direct.reduced.basis, padded.reduced.basis);
    assert_eq!(direct.reduced.targets, padded.reduced.targets);
    assert_eq!(direct.reduced.candidates, padded.reduced.candidates);
    assert_eq!(direct.reduced.matrix.len(), padded.reduced.matrix.len());
    for (direct, padded) in direct
        .reduced
        .matrix
        .iter()
        .flatten()
        .zip(padded.reduced.matrix.iter().flatten())
    {
        assert!((direct - padded).together().cancel().is_zero());
    }
}
