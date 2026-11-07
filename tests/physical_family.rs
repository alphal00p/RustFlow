use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::{reduction::Reduction, *};

fn tadpole(mass: Atom, scale: Atom) -> IntegralFamily {
    IntegralFamily {
        name: "physical_family_identity_test".into(),
        loops: vec!["k".into()],
        external: vec![],
        external_gram: vec![],
        propagators: vec![Propagator {
            constant: -&scale * mass,
            scalar_products: vec![scale],
        }],
        physical_propagators: 1,
        epsilon: symbol!("physical_family_test::epsilon"),
        dimension: 4,
    }
}
fn table(
    rules: BTreeMap<Integral, BTreeMap<Integral, Atom>>,
    conditions: Vec<Atom>,
) -> TableBackend {
    TableBackend {
        name: "physical_family_exact_tadpole".into(),
        reduction: Reduction {
            rules,
            residuals: vec![Integral(vec![1])],
            nonzero_conditions: conditions,
        },
    }
}
#[test]
fn independent_partial_derivatives_close_even_when_their_sum_cancels() {
    let (s, t) = symbol!("physical_family_test::s", "physical_family_test::t");
    let mass = Atom::var(s) - Atom::var(t);
    let family = tadpole(mass.clone(), Atom::one());
    let coefficient = (Atom::one() - Atom::var(family.epsilon)) / &mass;
    let backend = table(
        BTreeMap::from([(
            Integral(vec![2]),
            BTreeMap::from([(Integral(vec![1]), coefficient.clone())]),
        )]),
        vec![mass],
    );
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[s, t],
        &backend,
        &FlowOptions::default(),
        "one mass sheet",
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(prepared.basis(), &[Integral(vec![1])]);
    let derivatives = &prepared.flow().system().derivatives;
    assert!(
        (&derivatives[&s][0][0] - &coefficient)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        (&derivatives[&t][0][0] + &coefficient)
            .together()
            .cancel()
            .is_zero()
    );
    assert!(
        (&derivatives[&s][0][0] + &derivatives[&t][0][0])
            .together()
            .cancel()
            .is_zero()
    );
}
#[test]
fn hidden_normalization_parameters_are_rejected_before_reduction() {
    struct ForbiddenBackend;
    impl ReductionBackend for ForbiddenBackend {
        fn identity(&self) -> String {
            "must not run".into()
        }
        fn reduce(&self, _: &IntegralFamily, _: &[Integral], _: &RunContext) -> Result<Reduction> {
            panic!("hidden family parameter must fail before reduction")
        }
    }
    let (s, u) = symbol!("physical_family_test::s", "physical_family_test::hidden_u");
    let family = tadpole(Atom::one(), Atom::var(u));
    let result = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[s],
        &ForbiddenBackend,
        &FlowOptions::default(),
        "one sheet",
        &RunContext::default(),
    );
    assert!(matches!(result, Err(Error::InvalidInput(_))));
}
#[test]
fn identical_connections_with_distinct_family_normalizations_have_distinct_identities() {
    let variable = symbol!("physical_family_test::s");
    let backend = table(BTreeMap::new(), vec![]);
    let prepare = |scale| {
        PreparedPhysicalFamily::new(
            &tadpole(Atom::one(), Atom::num(scale)),
            &[Integral(vec![1])],
            &[variable],
            &backend,
            &FlowOptions::default(),
            "one sheet",
            &RunContext::default(),
        )
        .unwrap()
    };
    let first = prepare(1);
    let second = prepare(2);
    assert_eq!(
        first.flow().system().derivatives,
        second.flow().system().derivatives
    );
    assert_ne!(
        first.flow().identity().key(),
        second.flow().identity().key()
    );
}

#[test]
fn physical_seed_checkpoints_refined_samples_and_repairs_only_when_reuse_is_disabled() {
    let mass = symbol!("physical_family_test::checkpoint_mass");
    // The factor of two in the propagator must survive the unit projections.
    let family = tadpole(Atom::var(mass), Atom::num(2));
    let backend = RustRedBackend::default();
    let context = RunContext::default();
    let mut options = FlowOptions {
        digits: 12,
        ..Default::default()
    };
    let prepared = PreparedPhysicalFamily::new(
        &family,
        &[Integral(vec![1])],
        &[mass],
        &backend,
        &options,
        "positive tadpole mass",
        &context,
    )
    .unwrap();
    let point = BTreeMap::from([(mass, Atom::num(3))]);
    let extra_digits = 4;
    let seed = |cache: &mut RustFlowCache, options: &FlowOptions| {
        prepared.seed_cache(cache, &point, 0, extra_digits, options, &backend, &context)
    };
    let direct = seed(&mut RustFlowCache::default(), &options).unwrap();
    let directory = std::env::temp_dir().join(format!(
        "rustflow-physical-seed-samples-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&directory);
    options.sample_cache_directory = Some(directory.clone());
    let cached = seed(&mut RustFlowCache::default(), &options).unwrap();
    let checkpoints = std::fs::read_dir(&directory)
        .unwrap()
        .map(|entry| {
            let path = entry.unwrap().path();
            assert!(
                path.file_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .starts_with("sample-")
            );
            let bytes = std::fs::read(&path).unwrap();
            (path, bytes)
        })
        .collect::<BTreeMap<_, _>>();
    // Both independent epsilon grids must be checkpointed, including the
    // higher precision/order samples that establish the boundary's accuracy.
    let minimum_samples = (0..2)
        .map(|refinement| {
            epsilon::EpsilonFitProfile::new(
                options.digits + extra_digits,
                options.guard_digits,
                -2,
                0,
                refinement,
            )
            .unwrap()
            .samples
            .len()
        })
        .sum::<usize>();
    assert!(checkpoints.len() >= minimum_samples);
    let (used, pristine) = checkpoints.first_key_value().unwrap();
    let mut corrupt = pristine.clone();
    *corrupt.last_mut().unwrap() ^= 1;
    std::fs::write(used, &corrupt).unwrap();
    let mut refused = RustFlowCache::default();
    assert!(matches!(seed(&mut refused, &options), Err(Error::Cache(_))));
    assert!(refused.entries().is_empty());
    assert_eq!(std::fs::read(used).unwrap(), corrupt);

    options.reuse_samples = false;
    let repaired = seed(&mut RustFlowCache::default(), &options).unwrap();
    assert_eq!(std::fs::read(used).unwrap(), *pristine);
    options.reuse_samples = true;
    let warm = seed(&mut RustFlowCache::default(), &options).unwrap();
    assert_eq!(warm.coefficients, cached.coefficients);
    assert_eq!(
        warm.accuracy.comparison_errors(),
        cached.accuracy.comparison_errors()
    );
    assert_eq!(
        std::fs::read_dir(&directory).unwrap().count(),
        checkpoints.len()
    );
    for (path, bytes) in &checkpoints {
        assert_eq!(std::fs::read(path).unwrap(), *bytes);
    }

    let p = Precision::decimal(80).unwrap();
    for boundary in [&direct, &cached, &repaired, &warm] {
        boundary.validate().unwrap();
        assert_eq!(boundary.identity.key(), prepared.flow().identity().key());
        assert!(matches!(&boundary.point, transport_cache::CachedPoint::Exact(p) if p == &point));
        assert_eq!(boundary.kind, transport_cache::PointKind::Physical);
        assert_eq!(boundary.range, direct.range);
        assert_eq!(
            boundary.accuracy.verified_digits(),
            options.digits + extra_digits
        );
        assert_eq!(boundary.accuracy.provenance(), direct.accuracy.provenance());
        for (found, expected) in boundary.coefficients.iter().zip(&direct.coefficients) {
            assert_eq!(found.len(), prepared.basis().len());
            assert!(p.close(&found[0], &expected[0], options.digits));
        }
        // I_1 = -Gamma(epsilon-1) * 3^(1-epsilon) / 2 in the
        // ordinary d^D l/(i*pi^(D/2)) convention, without extra factors.
        assert!(p.close(&boundary.coefficients[0][0], &p.zero(), options.digits));
        assert!(p.close(
            &boundary.coefficients[1][0],
            &p.scale(&p.i(3), 1, 2),
            options.digits
        ));
        let gamma = ComplexFloat::new(p.real(0).euler(), p.real(0));
        let finite = p.scale(&p.sub(&p.sub(&p.i(1), &gamma), &p.log(&p.i(3))), 3, 2);
        assert!(p.close(&boundary.coefficients[2][0], &finite, options.digits));
    }
    std::fs::remove_dir_all(directory).unwrap();
}
