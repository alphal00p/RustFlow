//! Typed source-origin admission, separately from endpoint numerical accuracy.
use std::collections::BTreeMap;
use symbolica::prelude::*;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::guarded::{GuardedMeasureIdentity, IndexRole};
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
        measure: format!("actual assigned massless sunset input={input}; cut={cut:?}"),
        support: "positive separated chemical endpoints; source-origin admission tested separately"
            .into(),
        orientation: "future occupied momenta from exact family routing".into(),
        normalization: "unscaled source identities, original target phases retained".into(),
        branch: "formal physical factors; no origin-zero identity from this description alone"
            .into(),
        deformation: "all uncut physical quadratics D-eta; shells and occupations fixed".into(),
    }
}
fn source_origin<const N: usize>(cuts: &[usize]) {
    let input = input().prepare().unwrap();
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
    source_origin::<7>(&[0]);
}
#[test]
fn fully_compact_massless_origin_keeps_polynomial_and_storage_domains() {
    source_origin::<12>(&[0, 1]);
}
