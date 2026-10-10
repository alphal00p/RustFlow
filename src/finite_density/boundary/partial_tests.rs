use super::*;

#[test]
fn actual_partial_regions_use_bound_virtual_zeros_and_ordinary_hard_seeds() {
    use crate::finite_density::DensityInput;
    use crate::finite_density::partial_origin::continuation::PartialContinuation;
    let input: DensityInput = serde_json::from_str(include_str!(
        "../../../examples/finite_density/massless_three_loop_chain.json"
    ))
    .unwrap();
    let input = input.prepare().unwrap();
    let family = input
        .occupied_cut(&[0, 3], 4096)
        .unwrap()
        .at_physical_masses();
    let context = RunContext::default();
    let options = FlowOptions {
        digits: 20,
        guard_digits: 24,
        mass_mode: crate::MassMode::All,
        ..Default::default()
    };
    let backend = RustRedBackend {
        bubble_subloops: false,
        ..Default::default()
    };
    let epsilon_symbol = symbol!("rustflow_occupied::epsilon");
    let epsilon = Rational::from((-5, 4));
    let p = Precision::decimal(44).unwrap();
    let parameters = HashMap::from_iter([(Atom::var(epsilon_symbol), p.rational(&epsilon))]);
    let ordinary = family
        .region_family(epsilon_symbol, options.dimension)
        .unwrap();
    let compact = (0..family.loops())
        .map(|i| family.shells().iter().any(|s| s.loop_index == i))
        .collect::<Vec<_>>();
    let regions = crate::cut_regions::enumerate_compact(&ordinary, &compact, &context).unwrap();
    let master = family.targets()[0].keys().next().unwrap();
    let distributions = family
        .shells()
        .iter()
        .map(|s| OccupiedBoundaryDistribution {
            source_loop_index: s.loop_index,
            shell: CompactShell {
                mass_squared: Rational::zero(),
                chemical_potential: s.chemical_potential.clone(),
            },
            cut_index: master.0[s.physical_slot],
            upper_index: master.0[s.upper_slot],
            lower_index: master.0[s.lower_slot],
        })
        .collect::<Vec<_>>();
    let mut target = crate::Integral(master.0[..family.input_slots()].to_vec());
    for shell in family.shells() {
        target.0[shell.physical_slot] = 0;
    }
    for shifted in [vec![1, 2], vec![2, 4]] {
        let proof = PartialContinuation::new(
            &input,
            &family,
            &shifted,
            Default::default(),
            Default::default(),
            &context,
        )
        .unwrap();
        let integrated =
            IntegratedOccupiedBoundary::new(&backend, &options, &context, Default::default())
                .unwrap()
                .with_partial_flow_origin(&proof);
        let mask = (0..family.input_slots())
            .map(|i| shifted.contains(&i))
            .collect::<Vec<_>>();
        let mut zero_proofs = 0;
        let mut nonzero_hard = 0;
        for region in &regions {
            let transformed = ordinary.transform_loops(&region.transformation).unwrap();
            let expansion =
                crate::regions::expand_region(&ordinary, &target, &mask, region, 2).unwrap();
            for expression in &expansion.coefficients {
                let inherited = partial_region_conditions(
                    expression,
                    &expansion.coordinates,
                    &compact,
                    transformed.external.len(),
                    Default::default(),
                    &context,
                )
                .unwrap();
                let projected = crate::integrand::projected_factor_region(
                    &expression.together().cancel(),
                    &expansion.coordinates,
                    &transformed,
                    &region.hard,
                    10000,
                )
                .unwrap();
                // The public entry cannot silently omit pre-projection domains.
                assert!(matches!(
                    integrated.evaluate_projected_with_provenance(
                        &projected,
                        &distributions,
                        &epsilon,
                        &parameters,
                        p
                    ),
                    Err(Error::InvalidInput(_))
                ));
                let value = integrated
                    .evaluate_projected_with_inherited_conditions(
                        &projected,
                        &distributions,
                        &epsilon,
                        &parameters,
                        p,
                        Some(&inherited),
                    )
                    .unwrap();
                assert!(p.finite(&value.value));
                zero_proofs += value.virtual_soft_certificates.len();
                for certificate in &value.virtual_soft_certificates {
                    assert_eq!(certificate["parent_continuation"], proof.source_identity());
                }
                if !projected.hard.template.loops.is_empty() && value.value != p.zero() {
                    nonzero_hard += 1;
                }
            }
        }
        assert!(
            zero_proofs > 0,
            "virtual-soft region path was not exercised"
        );
        assert!(nonzero_hard > 0, "ordinary hard seeds were not exercised");
    }
}

#[test]
fn raw_parent_domain_precedes_region_cancellation_and_virtual_zero() {
    let coordinates = [
        "partial_g00",
        "partial_g01",
        "partial_g02",
        "partial_g11",
        "partial_g12",
        "partial_g22",
        "partial_e0",
        "partial_e1",
        "partial_e2",
    ]
    .iter()
    .map(|s| Atom::parse(s, "partial_raw_test", Default::default()).unwrap())
    .collect::<Vec<_>>();
    let epsilon = symbol!("partial_raw_test::ep");
    let ep = Atom::var(epsilon);
    let expression =
        (ep.clone().pow(2) - Atom::one()) / (ep.clone() - Atom::one()) / &coordinates[5];
    let conditions = partial_region_conditions(
        &expression,
        &coordinates,
        &[true, true, false],
        1,
        Default::default(),
        &RunContext::default(),
    )
    .unwrap();
    assert!(!conditions.is_empty());
    assert!(
        crate::physical_conditions::validate_conditions_at(
            &conditions,
            epsilon,
            &BTreeMap::from([(epsilon, Atom::num(1))])
        )
        .is_err()
    );
    crate::physical_conditions::validate_conditions_at(
        &conditions,
        epsilon,
        &BTreeMap::from([(epsilon, Atom::num(2))]),
    )
    .unwrap();

    // The denominator is compact-only, even though polynomial cancellation
    // would erase it and an independent virtual polynomial loop is present.
    let compact_pole = (coordinates[0].clone().pow(2) - Atom::one())
        / (&coordinates[0] - Atom::one())
        * &coordinates[5];
    assert!(matches!(
        partial_region_conditions(
            &compact_pole,
            &coordinates,
            &[true, true, false],
            1,
            Default::default(),
            &RunContext::default()
        ),
        Err(Error::Unsupported(_))
    ));

    let run = RunContext::default();
    run.cancellation.cancel();
    assert!(matches!(
        partial_region_conditions(
            &expression,
            &coordinates,
            &[true, true, false],
            1,
            Default::default(),
            &run
        ),
        Err(Error::Cancelled)
    ));
}
