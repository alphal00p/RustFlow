//! A genuine two-master bubble basis with an epsilon-singular physical connection.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{
    reduction::{LinearCombination, Reduction},
    transport_cache::*,
    *,
};

fn fixture() -> (IntegralFamily, Symbol, Reduction, Vec<Integral>) {
    let m = symbol!("physical_epsilon_shear::m");
    let eps = symbol!("physical_epsilon_shear::eps");
    let mass = Atom::var(m);
    let epsilon = Atom::var(eps);
    let family = IntegralFamily {
        name: "equal_mass_bubble_evanescent_tadpole_basis".into(),
        loops: vec!["k".into()],
        external: vec!["p".into()],
        external_gram: vec![vec![Atom::num(-1)]],
        propagators: vec![
            Propagator {
                constant: -&mass,
                scalar_products: vec![Atom::one(), Atom::new()],
            },
            Propagator {
                constant: -&mass - Atom::one(),
                scalar_products: vec![Atom::one(), Atom::num(2)],
            },
        ],
        physical_propagators: 2,
        epsilon: eps,
        dimension: 4,
    };
    let bubble = Integral(vec![1, 1]);
    let t3 = Integral(vec![3, 0]);
    let denominator = 4 * &mass + Atom::one();
    // These are physical identities, independently checked by native IBP in
    // the next test. The bubble relation follows by integrating the total
    // derivative of (x-1/2)[m+x(1-x)]^-eps on 0..1; the tadpole relation follows
    // from Gamma(a-D/2)/Gamma(a). Equal masses give I21=I12.
    let mixed = LinearCombination::from([
        (bubble.clone(), (Atom::one() - 2 * &epsilon) / &denominator),
        (t3.clone(), 2 * &mass / (&epsilon * &denominator)),
    ]);
    let fourth = LinearCombination::from([(t3.clone(), -(Atom::one() + &epsilon) / (3 * &mass))]);
    let table = Reduction {
        rules: BTreeMap::from([
            (Integral(vec![2, 1]), mixed.clone()),
            (Integral(vec![1, 2]), mixed),
            (Integral(vec![4, 0]), fourth),
        ]),
        residuals: vec![bubble.clone(), t3.clone()],
        ..Default::default()
    };
    (family, m, table, vec![bubble, t3, Integral(vec![2, 1])])
}

#[test]
fn supplied_evanescent_bubble_basis_has_exact_native_ibp_certificate() {
    let (family, _, table, targets) = fixture();
    let native = RustRedBackend {
        symmetry_rules: true,
        bubble_subloops: false,
        ..Default::default()
    };
    let mut all = targets;
    all.extend(table.rules.keys().cloned());
    all.sort();
    all.dedup();
    let reduction = native
        .reduce(&family, &all, &RunContext::default())
        .unwrap();
    for (target, rhs) in &table.rules {
        let mut difference = reduction.expand(target).unwrap();
        for (integral, coefficient) in rhs {
            for (master, weight) in reduction.expand(integral).unwrap() {
                *difference.entry(master).or_default() -= coefficient * weight;
            }
        }
        assert!(
            difference
                .values()
                .all(|coefficient| coefficient.together().cancel().is_zero()),
            "native IBP did not certify {target:?}: {difference:?}"
        );
    }
}

#[test]
fn amf_seeded_sheared_bubble_cache_restores_physical_targets_and_restarts() {
    let (family, mass, table, targets) = fixture();
    let context = RunContext::default();
    let options = FlowOptions::default();
    let supplied = ScopedTableBackend::from_table(
        "native-certified bubble and tadpole identities",
        &family,
        table,
    )
    .unwrap();
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &targets,
        &[mass],
        &supplied,
        &options,
        "positive mass with Euclidean external invariant",
        &context,
    )
    .unwrap();
    let retained = PreparedPhysicalFamily::new(
        &family,
        &targets[..2],
        &[mass],
        &supplied,
        &FlowOptions {
            skip_reduction: true,
            ..options.clone()
        },
        "positive mass with Euclidean external invariant",
        &context,
    )
    .unwrap();
    assert_eq!(prepared.basis(), &targets[..2]);
    assert_eq!(retained.basis(), prepared.basis());
    assert_eq!(
        retained.flow().system().derivatives,
        prepared.flow().system().derivatives
    );
    let original_key = prepared.flow().identity().key().to_owned();
    let original_basis = prepared.basis().to_vec();
    let original_targets = prepared.target_reductions().to_vec();
    let original_guards = prepared.nonzero_conditions().to_vec();
    let prepared = prepared.with_epsilon_shearing(&context).unwrap();
    assert_eq!(prepared.epsilon_shearing().unwrap().weights(), &[-1, 0]);
    assert_eq!(prepared.basis(), original_basis);
    assert_eq!(prepared.target_reductions(), original_targets);
    assert_eq!(prepared.nonzero_conditions(), original_guards);
    assert_ne!(prepared.flow().identity().key(), original_key);

    let source = BTreeMap::from([(mass, Atom::num(Rational::from((1, 2))))]);
    let destination = BTreeMap::from([(mass, Atom::num(Rational::from((2, 3))))]);
    let target_range = EpsilonRange::new(-2, 0).unwrap();
    let required = prepared
        .required_master_range(&destination, target_range)
        .unwrap();
    // Restoring B=epsilon^-1 J_B, and reducing I21, both require J through +1.
    assert_eq!(required, EpsilonRange::new(-2, 1).unwrap());
    let native = RustRedBackend::default();
    let mut cache = RustFlowCache::default();
    let seed = prepared
        .seed_cache(
            &mut cache,
            &source,
            required.last,
            16,
            &options,
            &native,
            &context,
        )
        .unwrap();
    assert_eq!(seed.range, required);
    assert!(seed.accuracy.verified_digits() >= 36);
    let short = CachedBoundary {
        identity: seed.identity.clone(),
        point: seed.point.clone(),
        kind: seed.kind,
        range: target_range,
        coefficients: seed.coefficients[..3].to_vec(),
        accuracy: BoundaryAccuracy::supplied(
            seed.accuracy.verified_digits(),
            seed.accuracy.working_bits(),
            seed.accuracy.comparison_errors()[..3].to_vec(),
            "deliberately truncated verified AMF seed for missing-order refusal",
        )
        .unwrap(),
    };
    assert!(matches!(
        prepared.project_targets(&short, target_range, 20),
        Err(Error::InvalidInput(_))
    ));
    let policy = ScaledDistance {
        scales: BTreeMap::new(),
        admissible: |_: &CachedBoundary, _: &CachedPoint| Ok(true),
    };
    let transport_options = FlowOptions {
        digits: 28,
        ..options.clone()
    };
    let mut short_cache = RustFlowCache::default();
    short_cache.insert(short).unwrap();
    assert!(
        prepared
            .flow()
            .evaluate_to(
                &mut short_cache,
                &destination,
                required,
                &transport_options,
                &context,
                &policy,
            )
            .is_err()
    );
    assert_eq!(short_cache.entries().len(), 1);
    let transported = prepared
        .flow()
        .evaluate_to(
            &mut cache,
            &destination,
            required,
            &transport_options,
            &context,
            &policy,
        )
        .unwrap();
    assert!(transported.transport.is_some());
    let projected = prepared
        .project_targets(&transported.boundary, target_range, 20)
        .unwrap();
    // Fresh destination AMF calculation uses the native reducer and independently
    // refines epsilon grids, precision and order; it consumes no cache values.
    let direct_options = FlowOptions {
        digits: 28,
        guard_digits: 52,
        series_order: 112,
        ..options.clone()
    };
    let direct = solve_integrals(
        &family,
        &targets,
        &KinematicPoint(BTreeMap::from([(
            Atom::var(mass),
            destination[&mass].clone(),
        )])),
        0,
        &direct_options,
        &native,
        &context,
    )
    .unwrap();
    let p = Precision::decimal(100).unwrap();
    for (found, direct) in projected.iter().zip(&direct) {
        assert!(direct.verified_digits.unwrap() >= 28);
        assert!(direct.refinements > 0 && direct.validation_samples > 0);
        for power in -2..=0 {
            assert!(p.close(
                &found.coefficients[&power],
                &direct.coefficients[&power],
                20
            ));
        }
    }
    let directory = std::env::temp_dir().join(format!(
        "rustflow-sheared-bubble-restart-{}",
        std::process::id()
    ));
    cache.save(&directory).unwrap();
    let mut loaded = RustFlowCache::load(&directory).unwrap();
    let repeat = prepared
        .flow()
        .evaluate_to(
            &mut loaded,
            &destination,
            required,
            &transport_options,
            &context,
            &policy,
        )
        .unwrap();
    assert!(repeat.transport.is_none());
    assert_eq!(
        repeat.boundary.coefficients,
        transported.boundary.coefficients
    );
    let restarted_targets = prepared
        .project_targets(&repeat.boundary, target_range, 20)
        .unwrap();
    for (old, restarted) in projected.iter().zip(restarted_targets) {
        assert_eq!(old.coefficients, restarted.coefficients);
        assert_eq!(old.absolute_errors, restarted.absolute_errors);
    }
    std::fs::remove_dir_all(directory).unwrap();
}
