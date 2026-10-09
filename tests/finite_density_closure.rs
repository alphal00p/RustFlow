//! Deliberately invoked bounded native closure diagnostics for physical slots.
//! No oracle values or finite-density numerical amplitudes enter this test.
use std::collections::BTreeMap;
use std::time::Instant;

use symbolica::prelude::*;
use symbolica_amflow::RunContext;
use symbolica_amflow::finite_density::DensityInput;
use symbolica_amflow::finite_density::guarded::{GuardedDiscoveryOptions, GuardedMeasureIdentity};
use symbolica_amflow::finite_density::preparation::{
    PreparedWeightedSources, WeightedSourcePolicy,
};
use symbolica_amflow::finite_density::reduction::{
    WeightedClosureOptions, WeightedClosureOutcome, prepare_weighted_system,
};

fn physical_sources<const N: usize>(
    cut: &[usize],
    shifted: &[usize],
) -> (DensityInput, PreparedWeightedSources<N>) {
    physical_sources_with_policy(cut, shifted, WeightedSourcePolicy::default())
}

fn physical_sources_with_policy<const N: usize>(
    cut: &[usize],
    shifted: &[usize],
    policy: WeightedSourcePolicy,
) -> (DensityInput, PreparedWeightedSources<N>) {
    let input: DensityInput = serde_json::from_str(include_str!(
        "../examples/finite_density/massive_two_loop_sunset.json"
    ))
    .unwrap();
    physical_sources_from_input_with_policy(input, cut, shifted, policy)
}

fn physical_sources_from_input<const N: usize>(
    input: DensityInput,
    cut: &[usize],
    shifted: &[usize],
) -> (DensityInput, PreparedWeightedSources<N>) {
    physical_sources_from_input_with_policy(input, cut, shifted, WeightedSourcePolicy::default())
}

fn physical_sources_from_input_with_policy<const N: usize>(
    input: DensityInput,
    cut: &[usize],
    shifted: &[usize],
    policy: WeightedSourcePolicy,
) -> (DensityInput, PreparedWeightedSources<N>) {
    let family = input
        .prepare()
        .unwrap()
        .occupied_cut(cut, 16)
        .unwrap()
        .at_physical_masses();
    let eta = symbol!("sunset_closure::eta");
    let epsilon = symbol!("sunset_closure::epsilon");
    let preparation = family
        .guarded_sources_with_policy::<N>(
            epsilon,
            4,
            eta,
            shifted,
            16,
            vec![],
            GuardedMeasureIdentity {
                measure: format!(
                    "sunset occupied {cut:?}; assigned physical factors={:?}; d=4-2epsilon",
                    family.factors()
                ),
                support: "mu=1; each occupied energy 0<E<1; formal continued uncut +i0; no endpoint-contour certificate asserted".into(),
                orientation: "future-oriented occupied loop momenta; native Minkowski quadratic forms".into(),
                normalization: "unscaled source identities; target Wick phases retained; no numerical normalization used".into(),
                branch: "generic complex eta and dimension away from retained nonzero conditions; completion powers<=0".into(),
                deformation: format!("native D-eta only in physical slots {shifted:?}; shells and occupation fixed"),
            },
            policy,
        )
        .unwrap();
    (input, preparation)
}

fn diagnostic<const N: usize>(cut: &[usize], shifted: &[usize]) {
    let policy = std::env::var("RUSTFLOW_WEIGHTED_SOURCE_POLICY")
        .map(|value| value.parse::<WeightedSourcePolicy>().unwrap())
        .unwrap_or_default();
    let (input, preparation) = physical_sources_with_policy::<N>(cut, shifted, policy);
    let root = std::env::var_os("RUSTFLOW_WEIGHTED_CLOSURE_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-weighted-closure"));
    let report = root.join(format!("sunset-N{N}"));
    let budget = |name: &str, default: usize| {
        std::env::var(name)
            .ok()
            .map(|text| text.parse().unwrap())
            .unwrap_or(default)
    };
    let options = WeightedClosureOptions {
        max_rounds: budget("RUSTFLOW_WEIGHTED_ROUNDS", 4),
        max_frontier: budget("RUSTFLOW_WEIGHTED_FRONTIER", 256),
        max_requested: budget("RUSTFLOW_WEIGHTED_REQUESTED", 2048),
        discovery: GuardedDiscoveryOptions {
            max_depth: budget("RUSTFLOW_WEIGHTED_DEPTH", 2).try_into().unwrap(),
            max_domains: budget("RUSTFLOW_WEIGHTED_DOMAINS", 128),
            sample_seed: 0,
        },
        checkpoints: Some(report.clone()),
        ..Default::default()
    };
    std::fs::create_dir_all(&report).unwrap();
    std::fs::write(
        report.join("input.json"),
        serde_json::to_vec_pretty(&serde_json::json!({
            "input": input, "cut_slots": cut, "shifted_slots": shifted,
            "source_policy": policy.as_str(),
            "source_count": preparation.context.sources().sources().len(),
            "max_rounds": options.max_rounds, "max_depth": options.discovery.max_depth,
            "max_domains": options.discovery.max_domains,
            "max_frontier": options.max_frontier, "max_requested": options.max_requested,
            "source_context": format!("{:?}", preparation.context.sources()),
            "numerical_prediction": false
        }))
        .unwrap(),
    )
    .unwrap();
    let started = Instant::now();
    eprintln!("sunset N={N}: native bounded discovery starting");
    let run = RunContext {
        progress: Some(std::sync::Arc::new(move |event| {
            eprintln!("sunset N={N}: {event:?}");
        })),
        ..Default::default()
    };
    let result = prepare_weighted_system(
        &preparation.context,
        &preparation.targets,
        &preparation.deformation,
        options,
        &run,
    )
    .unwrap();
    let elapsed = started.elapsed().as_secs_f64();
    let (status, diagnostics, detail) = match &result {
        WeightedClosureOutcome::Closed(closed) => {
            if let Some(system) = closed.differential_system() {
                system.validate().unwrap();
            }
            (
                "closed",
                &closed.diagnostics,
                serde_json::json!({
                    "basis": closed.reduced.basis.iter().map(|i| &i.0).collect::<Vec<_>>(),
                    "matrix": closed.reduced.matrix.iter().map(|r| r.iter().map(|a|a.to_canonical_string()).collect::<Vec<_>>()).collect::<Vec<_>>(),
                    "conditions": closed.reduced.nonzero_conditions.iter().map(|a|a.to_canonical_string()).collect::<Vec<_>>(),
                    "uncovered_discovery": format!("{:?}", closed.discovery_unresolved)
                }),
            )
        }
        WeightedClosureOutcome::Unresolved(unresolved) => (
            "unresolved",
            &unresolved.diagnostics,
            serde_json::json!({
                "reason": unresolved.reason,
                "frontier": unresolved.provisional_frontier.iter().map(|a|a.to_vec()).collect::<Vec<_>>(),
                "application_unresolved": format!("{:?}", unresolved.unresolved),
                "discovery_unresolved": format!("{:?}", unresolved.discovery_unresolved),
                "conditions": unresolved.nonzero_conditions.iter().map(|a|a.to_canonical_string()).collect::<Vec<_>>()
            }),
        ),
    };
    let summary = serde_json::json!({
        "status": status, "runtime_seconds": elapsed, "slots": N,
        "rounds": diagnostics.rounds, "requested": diagnostics.requested,
        "provisional_sizes": diagnostics.provisional_sizes,
        "native_rules": diagnostics.native_rules,
        "native_rule_applications": diagnostics.native_rule_applications,
        "detail": detail, "numerical_prediction": false
    });
    std::fs::write(
        report.join("result.json"),
        serde_json::to_vec_pretty(&summary).unwrap(),
    )
    .unwrap();
    eprintln!(
        "sunset N={N}: {status}, {elapsed:.3}s; report={}",
        report.display()
    );
}

#[test]
fn physical_sunset_weighted_contexts_construct() {
    let (_, single) =
        physical_sources_with_policy::<7>(&[0], &[1, 2], WeightedSourcePolicy::TangentsThenLorentz);
    let (_, double) =
        physical_sources_with_policy::<9>(&[0, 1], &[2], WeightedSourcePolicy::TangentsThenLorentz);
    assert_eq!(single.targets.len(), 2);
    assert_eq!(double.targets.len(), 2);
    for compact_loop in 0..2 {
        assert!(double.context.sources().sources().iter().any(|source| {
            source
                .id
                .starts_with(&format!("angular-tangent/{compact_loop}/"))
        }));
        assert!(double.context.sources().sources().iter().any(|source| {
            source
                .id
                .starts_with(&format!("shell-tangent/{compact_loop}/"))
        }));
    }
    assert!(
        single
            .context
            .sources()
            .sources()
            .iter()
            .any(|source| source.id.starts_with("angular-tangent/0/"))
    );
}

#[test]
fn source_policy_is_bound_and_legacy_preserves_all_occupation_faces() {
    use symbolica_amflow::finite_density::guarded::IndexBounds;
    let (_, legacy) = physical_sources::<9>(&[0, 1], &[2]);
    let (_, active) =
        physical_sources_with_policy::<9>(&[0, 1], &[2], WeightedSourcePolicy::ActiveLorentz);
    assert!(
        legacy
            .context
            .sources()
            .measure_id()
            .contains("source presentation=legacy-lorentz")
    );
    assert_ne!(
        legacy.context.sources().measure_id(),
        active.context.sources().measure_id()
    );
    assert!(
        legacy
            .context
            .sources()
            .sources()
            .iter()
            .all(|s| !s.id.contains("tangent"))
    );
    assert!(
        legacy
            .context
            .sources()
            .sources()
            .iter()
            .filter(|s| s.id.starts_with("lorentz/"))
            .all(|s| s.domain.bounds()[5..]
                .iter()
                .all(|b| *b == IndexBounds::fixed(0)
                    || *b == IndexBounds::new(Some(1), None).unwrap()))
    );
    assert!(
        active
            .context
            .sources()
            .sources()
            .iter()
            .filter(|s| s.id.starts_with("lorentz/"))
            .any(|s| s.domain.bounds()[5..]
                .iter()
                .any(|b| *b == IndexBounds::new(Some(0), None).unwrap()))
    );
    assert!("unrecognized".parse::<WeightedSourcePolicy>().is_err());
}

#[test]
fn positive_mass_lower_support_is_zero_but_massless_limit_is_not_inferred() {
    let (mut input, positive) = physical_sources::<7>(&[0], &[1, 2]);
    let lower_surface = [1, 1, 1, 0, 0, 0, 1];
    assert!(positive.context.sources().is_zero(&lower_surface));
    assert!(!positive.context.sources().is_zero(&[1, 1, 1, 0, 0, 1, 0]));
    assert!(
        positive
            .context
            .sources()
            .measure_id()
            .contains("certified real empty support")
    );
    let found = positive
        .context
        .discover(vec![], [], Default::default())
        .unwrap();
    let bytes = found.program.encode(Default::default()).unwrap();
    let restored = positive.context.decode(&bytes, Default::default()).unwrap();
    let reduced = restored.reduce(lower_surface, Default::default()).unwrap();
    assert!(reduced.terms.is_empty() && reduced.unresolved.is_empty());

    input.edges[0].mass_squared = "0".into();
    let (_, massless) = physical_sources_from_input::<7>(input, &[0], &[1, 2]);
    assert!(!massless.context.sources().is_zero(&lower_surface));
    assert!(massless.context.decode(&bytes, Default::default()).is_err());
}

#[test]
fn admitted_indices_precede_zero_and_constant_terminal_classification() {
    let (_, prepared) = physical_sources::<9>(&[0, 1], &[2]);
    for (indices, coefficient) in [
        ([1, 1, 0, 1, 0, 0, 0, 0, 0], Atom::one()),
        ([1, 1, 0, 1, 0, 0, 0, 0, 0], Atom::zero()),
        ([0, 1, 0, 0, 0, -1, 0, 0, 0], Atom::one()),
    ] {
        assert!(
            prepare_weighted_system(
                &prepared.context,
                &[BTreeMap::from([(indices, coefficient)])],
                &prepared.deformation,
                Default::default(),
                &RunContext::default(),
            )
            .is_err()
        );
    }
}

#[test]
fn raw_target_domains_survive_exact_zero_integrals_and_reject_unknown_coefficients() {
    let (_, prepared) = physical_sources::<9>(&[0, 1], &[2]);
    let epsilon = symbol!("sunset_closure::epsilon");
    let eps = Atom::var(epsilon);
    let removable = (eps.pow(2) - Atom::one()) / (&eps - Atom::one());
    assert!(
        (removable.together().cancel() - (&eps + Atom::one()))
            .expand()
            .is_zero()
    );
    let zero_integral = [0, 1, 0, 0, 0, 0, 0, 0, 0];
    let result = prepare_weighted_system(
        &prepared.context,
        &[BTreeMap::from([(zero_integral, removable)])],
        &prepared.deformation,
        Default::default(),
        &RunContext::default(),
    )
    .unwrap();
    let WeightedClosureOutcome::Closed(closed) = result else {
        panic!("cut-zero target did not close")
    };
    assert!(closed.reduced.basis.is_empty());
    assert!(closed.reduced.nonzero_conditions.iter().any(|condition| {
        (condition - (&eps - Atom::one())).expand().is_zero()
            || (condition + (&eps - Atom::one())).expand().is_zero()
    }));
    for invalid in [
        Atom::var(symbol!("undeclared_weighted_target_coefficient")),
        Atom::var(symbol!("rustflow_occupied_indices::a_0")),
        eps.pow(Atom::num((1, 2))),
    ] {
        assert!(
            prepare_weighted_system(
                &prepared.context,
                &[BTreeMap::from([(zero_integral, invalid)])],
                &prepared.deformation,
                Default::default(),
                &RunContext::default(),
            )
            .is_err()
        );
    }
}

#[test]
#[ignore = "bounded native search diagnostic; run explicitly and inspect closed/unresolved report"]
fn physical_sunset_single_cut_native_closure() {
    diagnostic::<7>(&[0], &[1, 2]);
}

#[test]
#[ignore = "bounded native search diagnostic; run explicitly and inspect closed/unresolved report"]
fn physical_sunset_double_cut_native_closure() {
    diagnostic::<9>(&[0, 1], &[2]);
}

#[test]
#[ignore = "one-domain native certification diagnostic; inspect retained unresolved evidence"]
fn native_single_domain_recenter_diagnostic() {
    use symbolica_amflow::finite_density::guarded::{IndexBounds, IndexDomain};
    let policy = std::env::var("RUSTFLOW_WEIGHTED_SOURCE_POLICY")
        .map(|value| value.parse::<WeightedSourcePolicy>().unwrap())
        .unwrap_or(WeightedSourcePolicy::TangentsThenLorentz);
    let (input, prepared) = physical_sources_with_policy::<9>(&[0, 1], &[2], policy);
    // One concrete domain retained by the bounded physical search, rather than
    // its growing derivative frontier. No expected failure is presumed here.
    let domain = IndexDomain::new([
        IndexBounds::fixed(1),
        IndexBounds::fixed(1),
        IndexBounds::fixed(1),
        IndexBounds::new(Some(-3), Some(-2)).unwrap(),
        IndexBounds::fixed(0),
        IndexBounds::fixed(0),
        IndexBounds::fixed(0),
        IndexBounds::new(Some(1), None).unwrap(),
        IndexBounds::fixed(0),
    ])
    .unwrap();
    let started = Instant::now();
    let discovery = prepared
        .context
        .discover(
            vec![domain.clone()],
            [],
            GuardedDiscoveryOptions {
                max_depth: 3,
                max_domains: 8192,
                sample_seed: 0,
            },
        )
        .unwrap();
    let root = std::env::var_os("RUSTFLOW_WEIGHTED_CLOSURE_REPORT")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("rustflow-weighted-closure"));
    let report = root.join("single-domain-recenter");
    std::fs::create_dir_all(&report).unwrap();
    let encoded = discovery.program.encode(Default::default()).unwrap();
    // Decode replays the original guarded source proof; saved programs are not
    // accepted merely because the external diagnostic says they are valid.
    let replayed = prepared
        .context
        .decode(&encoded, Default::default())
        .unwrap();
    let target = [1, 1, 1, -2, 0, 0, 0, 1, 0];
    let reduced = replayed.reduce(target, Default::default()).unwrap();
    std::fs::write(report.join("native-program.bin"), encoded).unwrap();
    let result = serde_json::json!({
        "input": input, "source_policy": policy.as_str(),
        "source_context": format!("{:?}", prepared.context.sources()),
        "domain": format!("{domain:?}"), "target": target,
        "max_depth": 3, "max_domains": 8192, "sample_seed": 0,
        "discovery_unresolved": format!("{:?}", discovery.unresolved),
        "application": format!("{reduced:?}"),
        "runtime_seconds": started.elapsed().as_secs_f64(),
        "numerical_prediction": false,
    });
    std::fs::write(
        report.join("result.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    eprintln!(
        "single-domain diagnostic: {} retained discovery gaps; report={}",
        discovery.unresolved.len(),
        report.display()
    );
}
