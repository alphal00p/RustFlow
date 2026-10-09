//! Typed source-origin admission, separately from endpoint numerical accuracy.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::guarded::{
    GuardedApplicationFailure, GuardedMeasureIdentity, IndexRole,
};
use symbolica_amflow::finite_density::massless_endpoint::MasslessFlowEvidence;
use symbolica_amflow::finite_density::preparation::{WeightedSourceOptions, WeightedSourcePolicy};
use symbolica_amflow::finite_density::reduction::{
    WeightedClosureOutcome, prepare_weighted_system,
};

fn input() -> DensityInput {
    let mut input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    for edge in &mut input.edges {
        edge.mass_squared = "0".into();
    }
    input
}
fn identity(input: &str, cut: &[usize]) -> GuardedMeasureIdentity {
    GuardedMeasureIdentity {
        measure: format!("actual assigned massless family input={input}; cut={cut:?}"),
        support: "positive separated chemical endpoints; source-origin admission tested separately"
            .into(),
        orientation: "future occupied momenta from exact family routing".into(),
        normalization: "unscaled source identities, original target phases retained".into(),
        branch: "formal physical factors; no origin-zero identity from this description alone"
            .into(),
        deformation: "all uncut physical quadratics D-eta; shells and occupations fixed".into(),
    }
}
fn source_origin<const N: usize>(definition: DensityInput, cuts: &[usize]) {
    let input = definition.prepare().unwrap();
    let family = input.occupied_cut(cuts, 16).unwrap().at_physical_masses();
    let shifted = (0..family.physical_slots())
        .filter(|&slot| family.roles()[slot] == IndexRole::Ordinary)
        .collect::<Vec<_>>();
    let options = WeightedSourceOptions::default();
    let evidence = MasslessFlowEvidence::new(&input, &family, &shifted, options).unwrap();
    let epsilon = symbol!("massless_sources::epsilon");
    let eta = symbol!("massless_sources::eta");
    let plain = family
        .guarded_sources_with_options::<N>(
            epsilon,
            4,
            eta,
            &shifted,
            16,
            vec![],
            identity(input.identity(), cuts),
            options,
        )
        .unwrap();
    let continued = family
        .guarded_sources_with_massless_origin::<N>(
            epsilon,
            4,
            eta,
            &shifted,
            16,
            vec![],
            identity(input.identity(), cuts),
            options,
            &evidence,
        )
        .unwrap();
    assert_eq!(continued.context.physical_arity(), family.factors().len());
    // The native identity stores GuardedMeasureIdentity as JSON. Compare its
    // escaped string payload, retaining quotes/backslashes in structural data.
    let encoded_identity = serde_json::to_string(&evidence.source_identity()).unwrap();
    assert!(
        continued
            .context
            .sources()
            .measure_id()
            .contains(&encoded_identity[1..encoded_identity.len() - 1])
    );
    assert!(
        continued
            .context
            .sources()
            .measure_id()
            .contains(evidence.origin_identity())
    );
    assert!(
        continued
            .context
            .sources()
            .measure_id()
            .contains("joint dimensional origin zero jets")
    );
    assert!(
        !continued
            .context
            .sources()
            .measure_id()
            .contains("certified real empty support")
    );
    let encoded = plain
        .context
        .discover(vec![], [], Default::default())
        .unwrap()
        .program
        .encode(Default::default())
        .unwrap();
    assert!(
        continued
            .context
            .decode(&encoded, Default::default())
            .is_err()
    );
    let program = continued
        .context
        .discover(vec![], [], Default::default())
        .unwrap()
        .program;
    let replayed = continued
        .context
        .decode(
            &program.encode(Default::default()).unwrap(),
            Default::default(),
        )
        .unwrap();
    let bulk: [i64; N] = std::array::from_fn(|slot| i64::from(slot < family.physical_slots()));
    assert!(!continued.context.sources().is_zero(&bulk));
    for shell in family.shells() {
        let mut upper = bulk;
        upper[shell.upper_slot] = 1;
        assert!(!continued.context.sources().is_zero(&upper));
        for cut_power in [1, 2, 4] {
            for lower_power in [1, 2, 3] {
                for upper_power in [0, 1, 2] {
                    let mut lower = bulk;
                    lower[shell.physical_slot] = cut_power;
                    lower[shell.lower_slot] = lower_power;
                    lower[shell.upper_slot] = upper_power;
                    assert!(!plain.context.sources().is_zero(&lower));
                    assert!(continued.context.sources().is_zero(&lower));
                    let reduced = replayed.reduce(lower, Default::default()).unwrap();
                    assert!(reduced.terms.is_empty() && reduced.unresolved.is_empty());
                }
            }
        }
        let mut lower = bulk;
        lower[shell.lower_slot] = 1;
        let outcome = prepare_weighted_system(
            &continued.context,
            &[BTreeMap::from([(lower, Atom::one())])],
            &continued.deformation,
            Default::default(),
            &Default::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Closed(closed) = outcome else {
            panic!("typed lower-contact zero did not close")
        };
        assert!(closed.reduced.basis.is_empty());
        // Invalid inverse completion labels cannot exploit the source zero.
        lower[family.physical_slots()] = 1;
        assert!(!continued.context.sources().is_zero(&lower));
        assert!(continued.deformation.derivative(lower).is_err());
    }
    if N > family.factors().len() {
        let mut escaped = bulk;
        escaped[family.shells()[0].lower_slot] = 1;
        escaped[family.factors().len()] = 1;
        assert!(replayed.reduce(escaped, Default::default()).is_err());
    }
    let different = WeightedSourceOptions {
        policy: WeightedSourcePolicy::ActiveLorentz,
        ..options
    };
    assert!(
        family
            .guarded_sources_with_massless_origin::<N>(
                epsilon,
                4,
                eta,
                &shifted,
                16,
                vec![],
                identity(input.identity(), cuts),
                different,
                &evidence
            )
            .is_err()
    );
    let inverse = WeightedSourceOptions {
        positive_compact_energy_powers: true,
        ..options
    };
    assert!(
        family
            .guarded_sources_with_massless_origin::<N>(
                epsilon,
                4,
                eta,
                &shifted,
                16,
                vec![],
                identity(input.identity(), cuts),
                inverse,
                &evidence
            )
            .is_err()
    );
    assert!(
        family
            .guarded_sources_with_massless_origin::<N>(
                epsilon,
                4,
                eta,
                &[],
                16,
                vec![],
                identity(input.identity(), cuts),
                options,
                &evidence
            )
            .is_err()
    );
}

#[test]
fn one_virtual_massless_origin_requires_its_sealed_prescription() {
    source_origin::<7>(input(), &[0]);
}
#[test]
fn fully_compact_massless_origin_keeps_polynomial_and_storage_domains() {
    source_origin::<12>(input(), &[0, 1]);
}

#[test]
fn four_loop_e7_sources_bind_each_sealed_endpoint_variant_and_zero_origin() {
    let definition: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/chain_of_three_parallel_pairs.json"
    ))
    .unwrap();
    for cuts in [vec![0], vec![4], vec![0, 4]] {
        source_origin::<24>(definition.clone(), &cuts);
    }
}

#[test]
fn sealed_free_virtual_zero_policy_is_opt_in_and_preserves_every_domain_bound() {
    let definition: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/chain_of_three_parallel_pairs.json"
    ))
    .unwrap();
    let input = definition.prepare().unwrap();
    let epsilon = symbol!("free_virtual_source_test::epsilon");
    let eta = symbol!("free_virtual_source_test::eta");
    for cuts in [vec![0], vec![4], vec![0, 4]] {
        let family = input.occupied_cut(&cuts, 16).unwrap().at_physical_masses();
        let shifted = (0..family.physical_slots())
            .filter(|s| !cuts.contains(s))
            .collect::<Vec<_>>();
        let off = WeightedSourceOptions::default();
        let on = WeightedSourceOptions {
            free_virtual_zero_sectors: true,
            ..off
        };
        let off_proof = MasslessFlowEvidence::new(&input, &family, &shifted, off).unwrap();
        let proof = MasslessFlowEvidence::new(&input, &family, &shifted, on).unwrap();
        assert!(
            matches!(family.guarded_sources_with_options::<24>(epsilon,4,eta,&shifted,16,vec![],identity(input.identity(),&cuts),on),
            Err(symbolica_amflow::Error::InvalidInput(ref message)) if message.contains("bound sealed"))
        );
        assert!(
            family
                .guarded_sources_with_massless_origin::<24>(
                    epsilon,
                    4,
                    eta,
                    &shifted,
                    16,
                    vec![],
                    identity(input.identity(), &cuts),
                    on,
                    &off_proof
                )
                .is_err()
        );
        let baseline = family
            .guarded_sources_with_massless_origin::<24>(
                epsilon,
                4,
                eta,
                &shifted,
                16,
                vec![],
                identity(input.identity(), &cuts),
                off,
                &off_proof,
            )
            .unwrap();
        let admitted = family
            .guarded_sources_with_massless_origin::<24>(
                epsilon,
                4,
                eta,
                &shifted,
                16,
                vec![],
                identity(input.identity(), &cuts),
                on,
                &proof,
            )
            .unwrap();
        let replay = admitted
            .context
            .discover(vec![], [], Default::default())
            .unwrap()
            .program;
        let encoded = replay.encode(Default::default()).unwrap();
        assert!(
            baseline
                .context
                .decode(&encoded, Default::default())
                .is_err()
        );
        let replay = admitted
            .context
            .decode(&encoded, Default::default())
            .unwrap();
        let bulk: [i64; 24] = std::array::from_fn(|slot| i64::from(slot < family.physical_slots()));
        assert!(!admitted.context.sources().is_zero(&bulk));
        for support in proof
            .free_virtual_zero_supports(&family, &shifted, 16)
            .unwrap()
        {
            let encoded_support = serde_json::to_string(&support.identity()).unwrap();
            assert!(
                admitted
                    .context
                    .sources()
                    .measure_id()
                    .contains(&encoded_support[1..encoded_support.len() - 1])
            );
            let mut polynomial = bulk;
            for &slot in support.forced_nonpositive_slots() {
                polynomial[slot] = 0;
            }
            assert!(!baseline.context.sources().is_zero(&polynomial));
            assert!(admitted.context.sources().is_zero(&polynomial));
            let reduced = replay.reduce(polynomial, Default::default()).unwrap();
            assert!(reduced.terms.is_empty() && reduced.unresolved.is_empty());
            let mut restored = polynomial;
            restored[support.forced_nonpositive_slots()[0]] = 1;
            assert!(!admitted.context.sources().is_zero(&restored));
            for (slot, value) in [
                (family.physical_slots(), 1),
                (family.shells()[0].upper_slot, -1),
                (family.factors().len(), 1),
            ] {
                let mut invalid = polynomial;
                invalid[slot] = value;
                assert!(!admitted.context.sources().is_zero(&invalid));
                if slot == family.factors().len() {
                    assert!(replay.reduce(invalid, Default::default()).is_err());
                } else {
                    // Native replay reports uncovered or role-invalid labels as
                    // explicit residuals; physical-domain rejection belongs to
                    // the weighted-system preparation layer below.
                    let residual = replay.reduce(invalid, Default::default()).unwrap();
                    assert!(residual.terms.is_empty());
                    assert_eq!(residual.unresolved.len(), 1);
                    assert_eq!(residual.unresolved[0].integral, invalid);
                    assert_eq!(residual.unresolved[0].coefficient, Atom::one());
                    assert_eq!(
                        residual.unresolved[0].reason,
                        if value < 0 {
                            GuardedApplicationFailure::InvalidOccupation { axis: slot }
                        } else {
                            GuardedApplicationFailure::NoApplicableRule
                        }
                    );
                }
                assert!(
                    prepare_weighted_system(
                        &admitted.context,
                        &[BTreeMap::from([(invalid, Atom::one())])],
                        &admitted.deformation,
                        Default::default(),
                        &Default::default(),
                    )
                    .is_err()
                );
            }
        }
        assert!(
            matches!(family.guarded_sources_with_massless_origin::<24>(epsilon,4,eta,&shifted,1,vec![],identity(input.identity(),&cuts),on,&proof),
            Err(symbolica_amflow::Error::Limit(ref message)) if message.contains("flat"))
        );
    }
}
