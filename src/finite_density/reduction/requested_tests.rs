//! Independent regressions for the optional requested-ray/exact-point schedule.
use super::*;
use crate::finite_density::guarded::{
    GuardedIdentity, GuardedIdentityTerm, GuardedMeasureIdentity,
};

fn options() -> WeightedClosureOptions {
    WeightedClosureOptions {
        active_target_closure: true,
        requested_index_rays: true,
        search_frontier_sectors: true,
        max_reused_rules: 8192,
        guard_refinement: GuardRefinementOptions {
            max_passes: 0,
            ..Default::default()
        },
        requested_ray_point: Some(RequestedRayPointOptions {
            max_domains_per_ray: 1,
            max_domains_per_point: 1,
        }),
        ..Default::default()
    }
}

fn same(a: &Atom, b: &Atom) -> bool {
    (a - b).together().cancel().is_zero()
}

fn coupled_context(with_point_pivot: bool) -> (GuardedContext<2>, FixedShellDeformation<2>) {
    let a = symbol!("requested_coupled_a");
    let b = symbol!("requested_coupled_b");
    let eta = symbol!("requested_coupled_eta");
    let parameter = symbol!("requested_coupled_parameter");
    let roles = [IndexRole::Ordinary; 2];
    let domain = IndexDomain::new([
        IndexBounds::new(Some(1), None).unwrap(),
        IndexBounds::new(Some(1), None).unwrap(),
    ])
    .unwrap();
    let rows = [
        ("original-coupled-pivot", Atom::var(a) - Atom::var(b)),
        ("original-point-pivot", Atom::var(eta)),
    ]
    .into_iter()
    .take(if with_point_pivot { 2 } else { 1 })
    .map(|(id, coefficient)| GuardedIdentity {
        id: id.into(),
        domain: domain.clone(),
        terms: vec![
            GuardedIdentityTerm {
                shift: [0, 0],
                coefficient,
            },
            GuardedIdentityTerm {
                shift: [-1, 0],
                coefficient: Atom::num(-1),
            },
        ],
        nonzero_conditions: vec![Atom::var(parameter)],
    })
    .collect();
    let context = GuardedContext::new(
        GuardedMeasureIdentity {
            measure: "synthetic original coupled-guard corpus".into(),
            support: "positive source indices".into(),
            orientation: "formal".into(),
            normalization: "formal; no period".into(),
            branch: "physical parameter nonzero".into(),
            deformation: "first ordinary index".into(),
        },
        roles,
        [a, b],
        vec![eta, parameter],
        rows,
    )
    .unwrap();
    let deformation = FixedShellDeformation::new(
        eta,
        AuxiliaryConvention::NativeMinusEta,
        roles,
        [true, false],
    )
    .unwrap();
    (context, deformation)
}

#[test]
fn requested_coupled_guard_needs_actual_point_proof_and_retains_parent_gap() {
    for available in [true, false] {
        let (context, deformation) = coupled_context(available);
        let requested = BTreeSet::from([[2, 2]]);
        let mut configured = options();
        configured.discovery.max_domains = 2;
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let mut conditions = BTreeMap::new();
        let (found, completed, deferred) = super::requested::prepare(
            &context,
            &requested,
            &deformation,
            &configured,
            &mut None,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(completed, requested);
        assert!(deferred.is_empty());
        assert_eq!(diagnostics.requested_ray_point_allocated_domains, 2);
        assert_eq!(diagnostics.requested_ray_point_fallbacks, 1);
        assert!(found.program.native().terminals().is_empty());
        assert!(
            found
                .unresolved
                .iter()
                .any(|gap| gap.reason == GuardedUnresolvedReason::ExceptionalCondition)
        );
        let replayed = context
            .decode(
                &found.program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        let result = replayed.reduce([2, 2], Default::default()).unwrap();
        if available {
            assert!(result.rule_applications > 0);
            assert!(result.unresolved.iter().all(|leaf| leaf.integral != [2, 2]
                && leaf.reason == GuardedApplicationFailure::NoApplicableRule));
            for symbol in [
                symbol!("requested_coupled_eta"),
                symbol!("requested_coupled_parameter"),
            ] {
                assert!(
                    result
                        .nonzero_conditions
                        .iter()
                        .any(|condition| condition.contains_symbol(symbol))
                );
                assert!(
                    conditions
                        .values()
                        .any(|condition| condition.contains_symbol(symbol))
                );
            }
        } else {
            assert!(result.unresolved.iter().any(|leaf| leaf.integral == [2, 2]
                && matches!(
                    leaf.reason,
                    GuardedApplicationFailure::ConditionVanished { .. }
                )));
        }
    }
}

#[test]
fn requested_budget_defers_whole_transactions_and_clears_stale_terminals() {
    let (context, deformation) = coupled_context(true);
    let requested = BTreeSet::from([[2, 2], [3, 3]]);
    for budget in [1, 2] {
        let mut configured = options();
        configured.discovery.max_domains = budget;
        // A previous provisional terminal must not suppress native discovery.
        let mut retained = Some(
            context
                .discover(vec![], [[2, 2]], Default::default())
                .unwrap()
                .program,
        );
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let (found, completed, deferred) = super::requested::prepare(
            &context,
            &requested,
            &deformation,
            &configured,
            &mut retained,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut BTreeMap::new(),
            &RunContext::default(),
        )
        .unwrap();
        assert!(found.program.native().terminals().is_empty());
        if budget == 1 {
            assert!(completed.is_empty());
            assert_eq!(deferred, requested);
            assert_eq!(diagnostics.requested_ray_point_allocated_domains, 0);
            assert_eq!(diagnostics.requested_ray_point_fallbacks, 0);
        } else {
            assert_eq!(completed, BTreeSet::from([[2, 2]]));
            assert_eq!(deferred, BTreeSet::from([[3, 3]]));
            assert_eq!(diagnostics.requested_ray_point_allocated_domains, 2);
            assert_eq!(diagnostics.requested_ray_point_fallbacks, 1);
        }
        assert_eq!(diagnostics.requested_ray_point_deferred, deferred.len());
    }
    let point = [2, 2];
    let mut covered = Some(
        context
            .discover(
                vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],
                [],
                Default::default(),
            )
            .unwrap()
            .program,
    );
    let mut configured = options();
    configured.discovery.max_domains = 1; // Insufficient for any new transaction.
    let mut diagnostics = WeightedClosureDiagnostics::default();
    let (_, completed, deferred) = super::requested::prepare(
        &context,
        &BTreeSet::from([point]),
        &deformation,
        &configured,
        &mut covered,
        &mut BTreeMap::new(),
        &mut diagnostics,
        &mut BTreeMap::new(),
        &RunContext::default(),
    )
    .unwrap();
    assert_eq!(completed, BTreeSet::from([point]));
    assert!(deferred.is_empty());
    assert_eq!(diagnostics.requested_ray_point_allocated_domains, 0);
    let mut configured = options();
    configured
        .requested_ray_point
        .as_mut()
        .unwrap()
        .max_domains_per_ray = usize::MAX;
    let mut diagnostics = WeightedClosureDiagnostics::default();
    let result = super::requested::prepare(
        &context,
        &requested,
        &deformation,
        &configured,
        &mut None,
        &mut BTreeMap::new(),
        &mut diagnostics,
        &mut BTreeMap::new(),
        &RunContext::default(),
    );
    assert!(matches!(result, Err(Error::Limit(_))), "{result:?}");
    assert_eq!(diagnostics.requested_ray_point_calls, 0);
}

#[test]
fn requested_already_covered_point_still_requires_original_source_binding() {
    let (context, deformation) = coupled_context(true);
    let (different, _) = coupled_context(false);
    let point = [2, 2];
    let mut retained = Some(
        context
            .discover(
                vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],
                [],
                Default::default(),
            )
            .unwrap()
            .program,
    );
    assert!(
        retained
            .as_ref()
            .unwrap()
            .reduce(point, Default::default())
            .unwrap()
            .rule_applications
            > 0
    );
    let result = super::requested::prepare(
        &different,
        &BTreeSet::from([point]),
        &deformation,
        &options(),
        &mut retained,
        &mut BTreeMap::new(),
        &mut Default::default(),
        &mut BTreeMap::new(),
        &RunContext::default(),
    );
    assert!(matches!(result, Err(Error::Reduction(_))), "{result:?}");
}

/// Independently reconstruct one complete combination using the final native
/// program. Only exact NoRule leaf coefficients may cancel within this sum.
fn replay_sum<const N: usize>(
    program: &GuardedReductionProgram<N>,
    input: &BTreeMap<[i64; N], Atom>,
) -> BTreeMap<[i64; N], Atom> {
    let mut terms = BTreeMap::<[i64; N], Atom>::new();
    let mut leaves = BTreeMap::<[i64; N], Atom>::new();
    for (label, coefficient) in input {
        let result = program.reduce(*label, Default::default()).unwrap();
        for (label, value) in result.terms {
            *terms.entry(label).or_default() += coefficient * value;
        }
        for residual in result.unresolved {
            assert_eq!(residual.reason, GuardedApplicationFailure::NoApplicableRule);
            *leaves.entry(residual.integral).or_default() += coefficient * residual.coefficient;
        }
    }
    assert!(leaves.values().all(|value| same(value, &Atom::zero())));
    terms.retain(|_, value| {
        *value = value.together().cancel();
        !value.is_zero()
    });
    terms
}

#[test]
fn requested_schedule_closes_compact_targets_and_replays_every_derivative() {
    let (context, deformation) = super::tests::compact_context();
    // Independent radial outputs; neither target weights nor basis are supplied
    // by a known period. The second output includes a genuine surface integral.
    let targets = [
        BTreeMap::from([([1, 0], Atom::one())]),
        BTreeMap::from([([2, 0], Atom::num(2)), ([1, 1], Atom::num(-1))]),
    ];
    let directory = std::env::temp_dir().join(format!(
        "rustflow-requested-ball-audit-{}",
        std::process::id()
    ));
    assert!(!directory.exists());
    let result = prepare_weighted_system(
        &context,
        &targets,
        &deformation,
        WeightedClosureOptions {
            checkpoints: Some(directory.clone()),
            ..options()
        },
        &RunContext::default(),
    )
    .unwrap();
    let WeightedClosureOutcome::Closed(closed) = result else {
        panic!("compact requested schedule did not close: {result:?}");
    };
    assert!(closed.diagnostics.requested_ray_point_calls > 0);
    assert!(closed.diagnostics.requested_ray_point_allocated_domains > 0);
    closed.differential_system().unwrap().validate().unwrap();
    let replayed = context
        .decode(
            &closed.program.encode(Default::default()).unwrap(),
            Default::default(),
        )
        .unwrap();
    for (target, expected) in targets.iter().zip(&closed.reduced.targets) {
        let actual = replay_sum(&replayed, target);
        let expected = expected
            .iter()
            .map(|(label, value)| {
                (
                    [i64::from(label.0[0]), i64::from(label.0[1])],
                    value.clone(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        assert_eq!(actual.len(), expected.len());
        for (label, value) in actual {
            assert!(same(&value, &expected[&label]));
        }
    }
    for (row, label) in closed.reduced.basis.iter().enumerate() {
        let label = [i64::from(label.0[0]), i64::from(label.0[1])];
        let actual = replay_sum(&replayed, &deformation.derivative(label).unwrap());
        for (column, basis) in closed.reduced.basis.iter().enumerate() {
            let basis = [i64::from(basis.0[0]), i64::from(basis.0[1])];
            let value = actual.get(&basis).cloned().unwrap_or_default();
            assert!(same(&value, &closed.reduced.matrix[row][column]));
        }
        assert!(actual.keys().all(|label| {
            closed
                .reduced
                .basis
                .iter()
                .any(|basis| [i64::from(basis.0[0]), i64::from(basis.0[1])] == *label)
        }));
    }
    assert!(!closed.reduced.nonzero_conditions.is_empty());
    let metadata: serde_json::Value = serde_json::from_slice(
        &std::fs::read(directory.join("round-000-provisional.json")).unwrap(),
    )
    .unwrap();
    let bytes = std::fs::read(directory.join("round-000-provisional.bin")).unwrap();
    assert_eq!(
        metadata["program_blake3"],
        blake3::hash(&bytes).to_hex().to_string()
    );
    assert_eq!(metadata["discovery_schedule"], "requested-ray-point-v1");
    let transaction_binding = &metadata["active_state"]["requested_discovery"];
    let transaction_bytes =
        std::fs::read(directory.join(transaction_binding["path"].as_str().unwrap())).unwrap();
    assert_eq!(
        transaction_binding["blake3"],
        blake3::hash(&transaction_bytes).to_hex().to_string()
    );
    let submitted = metadata["active_state"]["submitted_requests"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| [row[0].as_i64().unwrap(), row[1].as_i64().unwrap()])
        .collect::<BTreeSet<_>>();
    let mut expected = targets
        .iter()
        .flat_map(|target| target.keys().copied())
        .collect::<BTreeSet<_>>();
    for root in expected.clone() {
        expected.extend(deformation.derivative(root).unwrap().into_keys());
    }
    assert_eq!(
        submitted, expected,
        "round zero must include actual raw derivative labels"
    );
    let completed = metadata["active_state"]["completed_requests"]
        .as_array()
        .unwrap()
        .iter()
        .map(|row| [row[0].as_i64().unwrap(), row[1].as_i64().unwrap()])
        .collect::<BTreeSet<_>>();
    assert_eq!(completed, submitted);
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn requested_initial_derivatives_count_against_history_before_discovery() {
    let (context, deformation) = super::tests::compact_context();
    let result = prepare_weighted_system(
        &context,
        &[BTreeMap::from([([1, 0], Atom::one())])],
        &deformation,
        WeightedClosureOptions {
            max_requested: 1,
            ..options()
        },
        &RunContext::default(),
    )
    .unwrap();
    let WeightedClosureOutcome::Unresolved(failed) = result else {
        panic!("{result:?}")
    };
    assert!(failed.reason.contains("historical-request"));
    assert_eq!(failed.diagnostics.requested_ray_point_calls, 0);
    assert_eq!(failed.diagnostics.requested_ray_point_allocated_domains, 0);
}

#[test]
fn requested_cancelled_checkpoint_never_publishes_a_complete_marker() {
    let (context, _) = super::tests::compact_context();
    let program = context
        .discover(vec![], [], Default::default())
        .unwrap()
        .program;
    let directory = std::env::temp_dir().join(format!(
        "rustflow-requested-cancelled-checkpoint-{}",
        std::process::id(),
    ));
    assert!(!directory.exists());
    let configured = WeightedClosureOptions {
        checkpoints: Some(directory.clone()),
        ..options()
    };
    let run = RunContext::default();
    run.cancellation.cancel();
    assert!(matches!(
        checkpoint(
            &configured,
            0,
            true,
            &program,
            &BTreeSet::new(),
            &[],
            None,
            &run
        ),
        Err(Error::Cancelled)
    ));
    assert!(!directory.join("round-000-closed.json").exists());
    assert!(!directory.join("round-000-closed.bin").exists());
    std::fs::remove_dir_all(directory).unwrap();
}

#[test]
fn requested_actual_sum_cancellation_keeps_conditions_and_never_cancels_work_failure() {
    let eta = symbol!("requested_sum_eta");
    let x = symbol!("requested_sum_x");
    let roles = [IndexRole::Ordinary; 2];
    let context = GuardedContext::new(
        GuardedMeasureIdentity {
            measure: "formal two-moment source; no period supplied".into(),
            support: "exact original point identity".into(),
            orientation: "formal".into(),
            normalization: "formal".into(),
            branch: "x nonzero".into(),
            deformation: "ordinary factor zero".into(),
        },
        roles,
        [symbol!("requested_sum_a"), symbol!("requested_sum_b")],
        vec![eta, x],
        vec![GuardedIdentity {
            id: "original-equal-moments".into(),
            domain: IndexDomain::new([IndexBounds::fixed(1), IndexBounds::fixed(0)]).unwrap(),
            terms: vec![
                GuardedIdentityTerm {
                    shift: [0, 0],
                    coefficient: Atom::one(),
                },
                GuardedIdentityTerm {
                    shift: [-1, 1],
                    coefficient: Atom::num(-1),
                },
            ],
            nonzero_conditions: vec![Atom::var(x)],
        }],
    )
    .unwrap();
    let deformation = FixedShellDeformation::new(
        eta,
        AuxiliaryConvention::NativeMinusEta,
        roles,
        [true, false],
    )
    .unwrap();
    let targets = [BTreeMap::from([
        ([1, 0], Atom::one()),
        ([0, 1], Atom::num(-1)),
    ])];
    let result = prepare_weighted_system(
        &context,
        &targets,
        &deformation,
        options(),
        &RunContext::default(),
    )
    .unwrap();
    let WeightedClosureOutcome::Closed(closed) = result else {
        panic!("{result:?}")
    };
    assert!(closed.reduced.basis.is_empty());
    assert!(closed.reduced.targets[0].is_empty());
    assert!(closed.reduced.candidates.is_empty());
    assert!(closed.diagnostics.retired_unresolved > 0);
    assert!(
        closed
            .reduced
            .nonzero_conditions
            .iter()
            .any(|condition| condition.contains_symbol(x))
    );
    let replayed = context
        .decode(
            &closed.program.encode(Default::default()).unwrap(),
            Default::default(),
        )
        .unwrap();
    assert!(replay_sum(&replayed, &targets[0]).is_empty());
    let failed = prepare_weighted_system(
        &context,
        &targets,
        &deformation,
        WeightedClosureOptions {
            application: GuardedReductionLimits {
                max_rule_applications: 0,
                ..Default::default()
            },
            ..options()
        },
        &RunContext::default(),
    )
    .unwrap();
    let WeightedClosureOutcome::Unresolved(failed) = failed else {
        panic!("{failed:?}")
    };
    assert!(
        failed
            .unresolved
            .iter()
            .any(|residual| residual.reason == GuardedApplicationFailure::WorkLimit)
    );
}

#[test]
fn requested_schedule_is_default_off_and_rejects_conflicting_options() {
    assert!(
        WeightedClosureOptions::default()
            .requested_ray_point
            .is_none()
    );
    let (context, deformation) = super::tests::compact_context();
    let targets = [BTreeMap::from([([1, 0], Atom::one())])];
    let mut cases = Vec::new();
    for case in 0..10 {
        let mut configured = options();
        match case {
            0 => configured.active_target_closure = false,
            1 => configured.requested_index_rays = false,
            2 => configured.search_frontier_sectors = false,
            3 => configured.max_reused_rules = 0,
            4 => configured.max_domains_per_residual = 1,
            5 => configured.guard_refinement.max_passes = 1,
            6 => configured.prioritize_requested_indices = true,
            7 => configured.discovery.max_domains = 0,
            8 => {
                configured
                    .requested_ray_point
                    .as_mut()
                    .unwrap()
                    .max_domains_per_ray = 0
            }
            9 => {
                configured
                    .requested_ray_point
                    .as_mut()
                    .unwrap()
                    .max_domains_per_point = 0
            }
            _ => unreachable!(),
        }
        cases.push(configured);
    }
    for (case, configured) in cases.into_iter().enumerate() {
        assert!(
            matches!(
                prepare_weighted_system(
                    &context,
                    &targets,
                    &deformation,
                    configured,
                    &RunContext::default(),
                ),
                Err(Error::InvalidInput(_))
            ),
            "configuration case {case} silently accepted incompatible scheduling"
        );
    }
}

#[test]
fn requested_schedule_preserves_admission_before_any_native_work() {
    let eta = symbol!("requested_admission_eta");
    let roles = [
        IndexRole::Ordinary,
        IndexRole::Ordinary,
        IndexRole::Occupation,
        IndexRole::Ordinary,
    ];
    let bounds = [
        IndexBounds::unbounded(),
        IndexBounds::new(None, Some(0)).unwrap(),
        IndexBounds::new(Some(0), None).unwrap(),
        IndexBounds::fixed(0),
    ];
    let admitted = IndexDomain::new(bounds).unwrap();
    let context = GuardedContext::new_with_physical_arity(
        GuardedMeasureIdentity {
            measure: "formal admission regression, no physical period".into(),
            support: "polynomial completion and nonnegative occupation".into(),
            orientation: "formal".into(),
            normalization: "formal".into(),
            branch: "original source guard".into(),
            deformation: "ordinary first factor".into(),
        },
        roles,
        std::array::from_fn(|i| symbol!(format!("requested_admission_a{i}"))),
        vec![eta],
        vec![GuardedIdentity {
            id: "formal-original-row".into(),
            domain: admitted.clone(),
            terms: vec![GuardedIdentityTerm {
                shift: [0; 4],
                coefficient: Atom::var(eta),
            }],
            nonzero_conditions: vec![],
        }],
        3,
    )
    .unwrap();
    let deformation = FixedShellDeformation::new(
        eta,
        AuxiliaryConvention::NativeMinusEta,
        roles,
        [true, false, false, false],
    )
    .unwrap()
    .with_physical_arity(3)
    .unwrap()
    .with_admitted_domain(admitted)
    .unwrap();
    for label in [[1, 1, 0, 0], [1, 0, -1, 0], [1, 0, 0, 1]] {
        // Even a zero coefficient cannot erase invalid raw input admission.
        for coefficient in [Atom::zero(), Atom::one()] {
            let result = prepare_weighted_system(
                &context,
                &[BTreeMap::from([(label, coefficient)])],
                &deformation,
                options(),
                &RunContext::default(),
            );
            assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");
        }
    }
    let overflow = prepare_weighted_system(
        &context,
        &[BTreeMap::from([([i64::MAX, 0, 0, 0], Atom::one())])],
        &deformation,
        options(),
        &RunContext::default(),
    );
    assert!(matches!(overflow, Err(Error::Limit(_))), "{overflow:?}");
}
