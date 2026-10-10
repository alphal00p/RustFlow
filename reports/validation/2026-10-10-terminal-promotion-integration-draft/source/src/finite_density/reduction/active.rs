//! Reachability from actual target sums, with native discovery and replay.
use super::*;

#[derive(Default)]
struct Sum<const N: usize> {
    terms: BTreeMap<[i64; N], Atom>,
    leaves: BTreeMap<[i64; N], Atom>,
    failures: Vec<GuardedAtomUnresolved<N>>,
}

/// Only uncovered leaves can cancel within this one physical linear
/// combination. A native work limit, pole or invalid index remains a failure.
fn reduce_sum<const N: usize>(
    program: &GuardedReductionProgram<N>,
    input: &BTreeMap<[i64; N], Atom>,
    deformation: &FixedShellDeformation<N>,
    limits: GuardedReductionLimits,
    diagnostics: &mut WeightedClosureDiagnostics,
    conditions: &mut BTreeMap<String, Atom>,
    run: &RunContext,
) -> Result<Sum<N>> {
    let mut result = Sum::default();
    for (indices, weight) in input {
        deformation.validate_integral(indices)?;
        if weight.is_zero() {
            continue;
        }
        run.cancellation.check()?;
        let (reduced, memoized) = program.reduce_memoized(*indices, limits)?;
        let memo_usage = program.unit_reduction_memo_usage()?;
        if memoized {
            diagnostics.unit_memo_hits += 1;
            diagnostics.memoized_rule_applications += reduced.rule_applications;
        } else if memo_usage.is_some() {
            diagnostics.unit_memo_misses += 1;
        }
        if let Some(u) = memo_usage {
            diagnostics.unit_memo_peak_owned_bytes =
                diagnostics.unit_memo_peak_owned_bytes.max(u.owned_bytes);
            diagnostics.unit_memo_peak_entries = diagnostics.unit_memo_peak_entries.max(u.entries);
            diagnostics.unit_memo_peak_output_terms =
                diagnostics.unit_memo_peak_output_terms.max(u.output_terms);
            diagnostics.unit_memo_peak_conditions =
                diagnostics.unit_memo_peak_conditions.max(u.conditions);
            diagnostics.unit_memo_peak_polynomial_terms = diagnostics
                .unit_memo_peak_polynomial_terms
                .max(u.polynomial_terms);
        }
        diagnostics.native_rule_applications += reduced.rule_applications;
        retain_conditions(conditions, reduced.nonzero_conditions);
        for (label, coefficient) in reduced.terms {
            deformation.validate_integral(&label)?;
            *result.terms.entry(label).or_insert_with(Atom::new) += weight * coefficient;
        }
        for mut residual in reduced.unresolved {
            deformation.validate_integral(&residual.integral)?;
            residual.coefficient = (weight * residual.coefficient).together().cancel();
            if residual.reason == GuardedApplicationFailure::NoApplicableRule {
                *result
                    .leaves
                    .entry(residual.integral)
                    .or_insert_with(Atom::new) += residual.coefficient;
            } else {
                result.failures.push(residual);
            }
        }
    }
    for map in [&mut result.terms, &mut result.leaves] {
        map.retain(|_, coefficient| {
            *coefficient = coefficient.together().cancel();
            !coefficient.is_zero()
        });
    }
    Ok(result)
}

fn residuals<const N: usize>(sum: Sum<N>) -> Vec<GuardedAtomUnresolved<N>> {
    let mut result = sum.failures;
    result.extend(
        sum.leaves
            .into_iter()
            .map(|(integral, coefficient)| GuardedAtomUnresolved {
                integral,
                coefficient,
                reason: GuardedApplicationFailure::NoApplicableRule,
            }),
    );
    result
}

pub(super) fn prepare<const N: usize>(
    context: &GuardedContext<N>,
    targets: &[BTreeMap<[i64; N], Atom>],
    deformation: &FixedShellDeformation<N>,
    options: WeightedClosureOptions,
    mut conditions: BTreeMap<String, Atom>,
    run: &RunContext,
) -> Result<WeightedClosureOutcome<N>> {
    let roots = targets
        .iter()
        .flat_map(|target| target.iter())
        .filter(|(_, coefficient)| !coefficient.is_zero())
        .map(|(indices, _)| *indices)
        .collect::<BTreeSet<_>>();
    let mut requested = roots.clone();
    if options.requested_ray_point.is_some() {
        for root in &roots {
            requested.extend(deformation.derivative(*root)?.keys().copied());
        }
    }
    let mut history = BTreeSet::new();
    let mut diagnostics = WeightedClosureDiagnostics::default();
    let mut learned = Vec::new();
    let mut retained = None;
    let mut gaps = BTreeMap::new();
    let mut last_frontier = Vec::new();
    let mut last_residuals = Vec::new();
    let mut last_discovery = Vec::new();
    for round in 0..options.max_rounds {
        run.cancellation.check()?;
        if history.union(&requested).count() > options.max_requested {
            return Ok(unclosed(
                "active weighted historical-request budget exhausted",
                last_frontier,
                last_residuals,
                last_discovery,
                conditions,
                diagnostics,
            ));
        }
        let domains = discovery_domains(
            &requested,
            &deformation.roles,
            options.split_ordinary_zero_faces,
            options.requested_index_rays,
        )?
        .into_iter()
        .filter_map(|domain| domain.intersection(deformation.admitted_domain()))
        .collect::<Vec<_>>();
        // No earlier provisional stopping set is allowed to hide descendants.
        let (mut found, completed_requests, deferred_requests) =
            if options.requested_ray_point.is_some() {
                for label in &requested {
                    deformation.validate_integral(label)?;
                }
                discover_direct_zeros(
                    context,
                    &[],
                    &requested,
                    deformation,
                    &options,
                    &mut retained,
                    &mut gaps,
                    &mut diagnostics,
                    run,
                )?;
                requested::prepare(
                    context,
                    &requested,
                    deformation,
                    &options,
                    &mut retained,
                    &mut gaps,
                    &mut diagnostics,
                    &mut conditions,
                    run,
                )?
            } else {
                let found = discover_with_refinement(
                    context,
                    &domains,
                    &[],
                    &requested,
                    deformation,
                    &options,
                    &mut learned,
                    &mut retained,
                    &mut gaps,
                    &mut diagnostics,
                    &mut conditions,
                    run,
                )?;
                (found, requested.clone(), BTreeSet::new())
            };
        // New immutable epoch after discovery; no memo survives a previous round.
        found.program = found
            .program
            .with_unit_reduction_memo(options.unit_reduction_memo)?;
        history.extend(&completed_requests);
        diagnostics.rounds = round + 1;
        diagnostics.requested = requested.len();
        diagnostics.historical_requested = history.len();
        diagnostics.native_rules = found.program.native().rules().len();
        diagnostics.uncovered_discovery_domains = found.unresolved.len();
        let mut needed = roots.clone();
        needed.extend(&deferred_requests);
        let mut frontier = BTreeSet::new();
        let mut pending = BTreeSet::new();
        let mut failures = Vec::new();
        let mut uncovered = Vec::new();
        for target in targets {
            let sum = reduce_sum(
                &found.program,
                target,
                deformation,
                options.application,
                &mut diagnostics,
                &mut conditions,
                run,
            )?;
            pending.extend(sum.terms.keys().copied());
            pending.extend(sum.leaves.keys().copied());
            failures.extend(sum.failures);
            uncovered.extend(sum.leaves.into_iter().map(|(integral, coefficient)| {
                GuardedAtomUnresolved {
                    integral,
                    coefficient,
                    reason: GuardedApplicationFailure::NoApplicableRule,
                }
            }));
        }
        let mut deferred = !deferred_requests.is_empty();
        while failures.is_empty() {
            let Some(leaf) = pending.pop_first() else {
                break;
            };
            if !frontier.insert(leaf) {
                continue;
            }
            needed.insert(leaf);
            if frontier.len() > options.max_frontier {
                break;
            }
            // Every prospective terminal must receive a native discovery
            // attempt, including auxiliary-constant sectors.
            if !history.contains(&leaf) {
                deferred = true;
                continue;
            }
            let derivative = deformation.derivative(leaf)?;
            needed.extend(derivative.keys().copied());
            if derivative.keys().any(|index| !history.contains(index)) {
                deferred = true;
                continue;
            }
            let sum = reduce_sum(
                &found.program,
                &derivative,
                deformation,
                options.application,
                &mut diagnostics,
                &mut conditions,
                run,
            )?;
            pending.extend(
                sum.terms
                    .keys()
                    .filter(|index| !frontier.contains(*index))
                    .copied(),
            );
            pending.extend(
                sum.leaves
                    .keys()
                    .filter(|index| !frontier.contains(*index))
                    .copied(),
            );
            failures.extend(sum.failures);
            uncovered.extend(sum.leaves.into_iter().map(|(integral, coefficient)| {
                GuardedAtomUnresolved {
                    integral,
                    coefficient,
                    reason: GuardedApplicationFailure::NoApplicableRule,
                }
            }));
        }
        let frontier = frontier.into_iter().collect::<Vec<_>>();
        diagnostics.provisional_sizes.push(frontier.len());
        diagnostics.native_frontier_requests +=
            frontier.iter().filter(|i| !history.contains(*i)).count();
        let requested_discovery = if options.requested_ray_point.is_some() {
            options
                .checkpoints
                .as_ref()
                .map(|directory| {
                    let name = format!(
                        "requested-discovery-{:04}.json",
                        diagnostics.requested_ray_point_calls
                    );
                    let bytes = std::fs::read(directory.join(&name))?;
                    Ok::<_, Error>(serde_json::json!({"path":name,
                    "blake3":blake3::hash(&bytes).to_hex().to_string()}))
                })
                .transpose()?
        } else {
            None
        };
        let mut active_state = serde_json::json!({
            "requested_discovery":requested_discovery,
            "historical_requests":history.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
            "submitted_requests":requested.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
            "completed_requests":completed_requests.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
            "deferred_requests":deferred_requests.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
            "next_needed":needed.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
            "history_is_coverage_certificate":false,
        });
        run.cancellation.check()?;
        checkpoint(
            &options,
            round,
            false,
            &found.program,
            &needed,
            &frontier,
            Some(&active_state),
            run,
        )?;
        if let Some(directory) = &options.checkpoints {
            let path = directory.join("active-request-history.json");
            let temporary = directory.join("active-request-history.json.part");
            std::fs::write(&temporary, serde_json::to_vec_pretty(&serde_json::json!({
                "schema":1, "round":round, "historical_requests":history.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
                "current_requests":needed.iter().map(|i| i.to_vec()).collect::<Vec<_>>(),
                "historical_requests_are_additional_targets":false,
            })).map_err(|e| Error::Cache(e.to_string()))?)?;
            std::fs::rename(temporary, path)?;
        }
        if !failures.is_empty() {
            return Ok(unclosed(
                "native active weighted application failed",
                frontier,
                failures,
                found.unresolved,
                conditions,
                diagnostics,
            ));
        }
        if frontier.len() > options.max_frontier {
            return Ok(unclosed(
                "active weighted provisional-frontier budget exhausted",
                frontier,
                uncovered,
                found.unresolved,
                conditions,
                diagnostics,
            ));
        }
        run.emit(Progress::DifferentialClosure {
            round,
            requested: needed.len(),
            basis_size: frontier.len(),
            new_derivatives: needed.difference(&history).count(),
        })?;
        if !deferred {
            run.emit(Progress::Stage {
                name: "active weighted final complete source replay".into(),
            })?;
            run.cancellation.check()?;
            let program = if options.unit_reduction_memo.is_some() {
                found.program.with_promoted_terminals_replayed(
                    frontier.iter().copied(),
                    options.max_reused_rules,
                )?
            } else {
                found
                    .program
                    .with_terminals_replayed(frontier.iter().copied(), options.max_reused_rules)?
            };
            // Recompute every sum with this single final replayed program.
            let mut failed = Vec::new();
            let mut target_weights = Vec::new();
            run.emit(Progress::Stage {
                name: "active weighted final original target audit".into(),
            })?;
            for target in targets {
                let sum = reduce_sum(
                    &program,
                    target,
                    deformation,
                    options.application,
                    &mut diagnostics,
                    &mut conditions,
                    run,
                )?;
                target_weights.push(linear_combination(&sum.terms, deformation.physical_arity)?);
                failed.extend(residuals(sum));
            }
            let mut matrix = Vec::new();
            run.emit(Progress::Stage {
                name: "active weighted final basis derivative audit".into(),
            })?;
            for leaf in &frontier {
                let derivative = deformation.derivative(*leaf)?;
                let sum = reduce_sum(
                    &program,
                    &derivative,
                    deformation,
                    options.application,
                    &mut diagnostics,
                    &mut conditions,
                    run,
                )?;
                let mut row = vec![Atom::new(); frontier.len()];
                for (index, coefficient) in &sum.terms {
                    let column = frontier.binary_search(index).map_err(|_| {
                        Error::IncompleteReduction(
                            "active weighted derivative escaped its audited basis".into(),
                        )
                    })?;
                    row[column] = coefficient.clone();
                }
                matrix.push(row);
                failed.extend(residuals(sum));
            }
            if !failed.is_empty() {
                return Ok(unclosed(
                    "final active weighted target/derivative audit failed",
                    frontier,
                    failed,
                    found.unresolved,
                    conditions,
                    diagnostics,
                ));
            }
            let mut candidates = BTreeMap::new();
            let mut retired = Vec::new();
            // The deterministic ordered prefix is optional. It neither changes
            // the audited connection nor installs maps for unvisited history.
            let selected = options
                .max_history_candidate_maps
                .unwrap_or(usize::MAX)
                .min(history.len());
            diagnostics.history_candidate_selected = selected;
            diagnostics.history_candidate_unassessed = history.len() - selected;
            run.emit(Progress::Stage {
                name: format!(
                    "active weighted optional historical maps: {selected} selected, {} unassessed",
                    history.len() - selected
                ),
            })?;
            for label in history.iter().take(selected) {
                let reduced = reduce_sum(
                    &program,
                    &BTreeMap::from([(*label, Atom::num(1))]),
                    deformation,
                    options.application,
                    &mut diagnostics,
                    &mut conditions,
                    run,
                )?;
                diagnostics.history_candidate_assessed += 1;
                if reduced.failures.is_empty() && reduced.leaves.is_empty() {
                    diagnostics.history_candidate_maps += 1;
                    candidates.insert(
                        native_integral(label, deformation.physical_arity)?,
                        linear_combination(&reduced.terms, deformation.physical_arity)?,
                    );
                } else {
                    diagnostics.retired_unresolved += 1;
                    retired.push(serde_json::json!({"requested":label.to_vec(),
                        "unresolved":residuals(reduced).iter().map(|r| serde_json::json!({
                            "integral":r.integral.to_vec(), "coefficient":r.coefficient.to_canonical_string(),
                            "reason":format!("{:?}",r.reason)})).collect::<Vec<_>>()}));
                }
            }
            active_state["historical_candidate_collection"] = serde_json::json!({
                "max_history_candidate_maps": options.max_history_candidate_maps,
                "selection": "ascending stored-index lexicographic prefix",
                "selected": history.iter().take(selected).map(|i| i.to_vec()).collect::<Vec<_>>(),
                "unassessed": history.iter().skip(selected).map(|i| i.to_vec()).collect::<Vec<_>>(),
                "assessed_count": diagnostics.history_candidate_assessed,
                "certified_map_count": diagnostics.history_candidate_maps,
                "assessed_unresolved_count": diagnostics.retired_unresolved,
                "unassessed_are_unresolved": false,
                "history_is_coverage_certificate": false,
            });
            if let Some(directory) = &options.checkpoints {
                let path = directory.join("retired-unresolved.json");
                let temporary = directory.join("retired-unresolved.json.part");
                std::fs::write(
                    &temporary,
                    serde_json::to_vec_pretty(&retired).map_err(|e| Error::Cache(e.to_string()))?,
                )?;
                std::fs::rename(temporary, path)?;
            }
            run.emit(Progress::Stage {
                name: "active weighted final closed checkpoint".into(),
            })?;
            run.cancellation.check()?;
            checkpoint(
                &options,
                round,
                true,
                &program,
                &needed,
                &frontier,
                Some(&active_state),
                run,
            )?;
            let program = program.with_unit_reduction_memo(None)?;
            return Ok(WeightedClosureOutcome::Closed(WeightedReducedSystem {
                reduced: ReducedSystem {
                    basis: frontier
                        .iter()
                        .map(|i| native_integral(i, deformation.physical_arity))
                        .collect::<Result<_>>()?,
                    matrix,
                    targets: target_weights,
                    nonzero_conditions: conditions.into_values().collect(),
                    candidates,
                    transformations: Vec::new(),
                },
                variable: deformation.variable,
                roles: deformation.roles,
                program,
                discovery_unresolved: found.unresolved,
                diagnostics,
            }));
        }
        requested = needed;
        last_frontier = frontier;
        last_residuals = uncovered;
        last_discovery = found.unresolved;
        retained = Some(found.program.with_unit_reduction_memo(None)?);
    }
    Ok(unclosed(
        "active weighted closure round budget exhausted",
        last_frontier,
        last_residuals,
        last_discovery,
        conditions,
        diagnostics,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> WeightedClosureOptions {
        WeightedClosureOptions {
            active_target_closure: true,
            requested_index_rays: true,
            max_reused_rules: 8192,
            max_domains_per_residual: 32,
            ..Default::default()
        }
    }

    #[test]
    fn optional_unit_memo_preserves_full_weighted_final_audit() {
        let (context, deformation) = super::super::tests::compact_context();
        let targets = [BTreeMap::from([([1, 0], Atom::num(1))])];
        let make = |memo| {
            let outcome = prepare_weighted_system(
                &context,
                &targets,
                &deformation,
                WeightedClosureOptions {
                    unit_reduction_memo: memo,
                    max_history_candidate_maps: Some(0),
                    ..options()
                },
                &RunContext::default(),
            )
            .unwrap();
            let WeightedClosureOutcome::Closed(closed) = outcome else {
                panic!("{outcome:?}")
            };
            closed
        };
        let baseline = make(None);
        let limits = super::super::super::guarded::GuardedUnitMemoLimits::default();
        let memoized = make(Some(limits));
        assert_eq!(baseline.reduced.basis, memoized.reduced.basis);
        assert_eq!(baseline.reduced.matrix, memoized.reduced.matrix);
        assert_eq!(baseline.reduced.targets, memoized.reduced.targets);
        assert_eq!(
            baseline.reduced.nonzero_conditions,
            memoized.reduced.nonzero_conditions
        );
        assert_eq!(
            baseline.program.encode(Default::default()).unwrap(),
            memoized.program.encode(Default::default()).unwrap()
        );
        assert_eq!(
            baseline.diagnostics.native_rule_applications,
            memoized.diagnostics.native_rule_applications
        );
        assert!(memoized.diagnostics.unit_memo_hits > 0);
        assert!(memoized.diagnostics.unit_memo_peak_owned_bytes <= limits.max_owned_bytes);
        assert!(memoized.diagnostics.unit_memo_peak_conditions <= limits.max_conditions);
        assert!(
            memoized
                .program
                .unit_reduction_memo_usage()
                .unwrap()
                .is_none()
        );
    }

    #[test]
    fn unit_memo_requires_active_closure() {
        let (context, deformation) = super::super::tests::compact_context();
        let targets = [BTreeMap::from([([1, 0], Atom::num(1))])];
        let result = prepare_weighted_system(
            &context,
            &targets,
            &deformation,
            WeightedClosureOptions {
                unit_reduction_memo: Some(Default::default()),
                ..Default::default()
            },
            &RunContext::default(),
        );
        assert!(
            matches!(result, Err(Error::InvalidInput(message)) if message.contains("requires active closure"))
        );
    }

    #[test]
    fn bounded_history_maps_leave_the_required_connection_and_program_unchanged() {
        let (context, deformation) = super::super::tests::compact_context();
        let targets = [BTreeMap::from([([1, 0], Atom::num(1))])];
        let make = |limit| {
            let outcome = prepare_weighted_system(
                &context,
                &targets,
                &deformation,
                WeightedClosureOptions {
                    max_history_candidate_maps: limit,
                    ..options()
                },
                &RunContext::default(),
            )
            .unwrap();
            let WeightedClosureOutcome::Closed(closed) = outcome else {
                panic!("{outcome:?}")
            };
            closed
        };
        let all = make(None);
        assert!(all.diagnostics.historical_requested > 1);
        assert_eq!(all.diagnostics.history_candidate_unassessed, 0);
        for limit in [0, 1, usize::MAX] {
            let bounded = make(Some(limit));
            assert_eq!(bounded.reduced.basis, all.reduced.basis);
            assert_eq!(bounded.reduced.matrix, all.reduced.matrix);
            assert_eq!(bounded.reduced.targets, all.reduced.targets);
            assert_eq!(
                bounded.program.encode(Default::default()).unwrap(),
                all.program.encode(Default::default()).unwrap()
            );
            let diagnostics = &bounded.diagnostics;
            let selected = limit.min(diagnostics.historical_requested);
            assert_eq!(diagnostics.history_candidate_selected, selected);
            assert_eq!(diagnostics.history_candidate_assessed, selected);
            assert_eq!(
                diagnostics.history_candidate_maps,
                bounded.reduced.candidates.len()
            );
            assert_eq!(
                diagnostics.history_candidate_maps + diagnostics.retired_unresolved,
                selected
            );
            assert_eq!(
                selected + diagnostics.history_candidate_unassessed,
                diagnostics.historical_requested
            );
            for (label, map) in &bounded.reduced.candidates {
                assert_eq!(all.reduced.candidates.get(label), Some(map));
            }
            for condition in &bounded.reduced.nonzero_conditions {
                assert!(all.reduced.nonzero_conditions.contains(condition));
            }
        }
    }

    #[test]
    fn bounded_history_option_rejects_inactive_closure() {
        let (context, deformation) = super::super::tests::compact_context();
        let result = prepare_weighted_system(
            &context,
            &[BTreeMap::from([([1, 0], Atom::num(1))])],
            &deformation,
            WeightedClosureOptions {
                max_history_candidate_maps: Some(0),
                ..Default::default()
            },
            &RunContext::default(),
        );
        assert!(
            matches!(result, Err(Error::InvalidInput(message)) if message.contains("requires active"))
        );
    }

    #[test]
    fn skipped_history_does_not_suppress_required_work_failure() {
        let (context, deformation) = super::super::tests::compact_context();
        let result = prepare_weighted_system(
            &context,
            &[BTreeMap::from([([1, 0], Atom::num(1))])],
            &deformation,
            WeightedClosureOptions {
                max_history_candidate_maps: Some(0),
                application: GuardedReductionLimits {
                    max_rule_applications: 0,
                    ..Default::default()
                },
                ..options()
            },
            &RunContext::default(),
        )
        .unwrap();
        assert!(matches!(result, WeightedClosureOutcome::Unresolved(_)));
    }

    #[test]
    fn memo_cancellation_before_final_replay_cannot_publish_a_closed_result() {
        use std::sync::{Arc, Mutex};
        let (context, deformation) = super::super::tests::compact_context();
        let token = crate::CancellationToken::default();
        let captured = token.clone();
        let stages = Arc::new(Mutex::new(Vec::new()));
        let logged = stages.clone();
        let run = RunContext {
            cancellation: token,
            progress: Some(Arc::new(move |event| {
                if let Progress::Stage { name } = event {
                    logged.lock().unwrap().push(name.clone());
                    if name == "active weighted final complete source replay" {
                        captured.cancel();
                    }
                }
            })),
        };
        let result = prepare_weighted_system(
            &context,
            &[BTreeMap::from([([1, 0], Atom::num(1))])],
            &deformation,
            WeightedClosureOptions {
                max_history_candidate_maps: Some(0),
                unit_reduction_memo: Some(Default::default()),
                ..options()
            },
            &run,
        );
        assert!(matches!(result, Err(Error::Cancelled)));
        assert_eq!(
            *stages.lock().unwrap(),
            vec!["active weighted final complete source replay"]
        );
    }

    #[test]
    fn cancellation_before_final_replay_cannot_publish_a_closed_result() {
        use std::sync::{Arc, Mutex};
        let (context, deformation) = super::super::tests::compact_context();
        let token = crate::CancellationToken::default();
        let captured = token.clone();
        let stages = Arc::new(Mutex::new(Vec::new()));
        let logged = stages.clone();
        let run = RunContext {
            cancellation: token,
            progress: Some(Arc::new(move |event| {
                if let Progress::Stage { name } = event {
                    logged.lock().unwrap().push(name.clone());
                    if name == "active weighted final complete source replay" {
                        captured.cancel();
                    }
                }
            })),
        };
        let result = prepare_weighted_system(
            &context,
            &[BTreeMap::from([([1, 0], Atom::num(1))])],
            &deformation,
            WeightedClosureOptions {
                max_history_candidate_maps: Some(0),
                ..options()
            },
            &run,
        );
        assert!(matches!(result, Err(Error::Cancelled)));
        assert_eq!(
            *stages.lock().unwrap(),
            vec!["active weighted final complete source replay"]
        );
    }

    #[test]
    fn actual_sum_cancels_uncovered_leaves_but_never_work_failures() {
        let (context, deformation) = super::super::tests::compact_context();
        let program = context
            .discover(
                vec![IndexDomain::new([IndexBounds::fixed(2), IndexBounds::fixed(0)]).unwrap()],
                [],
                GuardedDiscoveryOptions {
                    max_depth: 3,
                    ..Default::default()
                },
            )
            .unwrap()
            .program;
        let first = program.reduce([2, 0], Default::default()).unwrap();
        assert!(!first.unresolved.is_empty());
        assert!(first.terms.is_empty());
        assert!(
            first
                .unresolved
                .iter()
                .all(|r| r.reason == GuardedApplicationFailure::NoApplicableRule)
        );
        let mut equation = BTreeMap::from([([2, 0], Atom::num(1))]);
        for residual in first.unresolved {
            assert_ne!(residual.integral, [2, 0]);
            *equation.entry(residual.integral).or_insert_with(Atom::new) -= residual.coefficient;
        }
        let mut conditions = BTreeMap::new();
        let result = reduce_sum(
            &program,
            &equation,
            &deformation,
            Default::default(),
            &mut Default::default(),
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert!(result.terms.is_empty() && result.leaves.is_empty() && result.failures.is_empty());
        assert!(!conditions.is_empty());
        let failed = reduce_sum(
            &program,
            &equation,
            &deformation,
            GuardedReductionLimits {
                max_rule_applications: 0,
                ..Default::default()
            },
            &mut Default::default(),
            &mut BTreeMap::new(),
            &RunContext::default(),
        )
        .unwrap();
        assert!(!failed.failures.is_empty());
        assert!(
            failed
                .failures
                .iter()
                .all(|r| r.reason == GuardedApplicationFailure::WorkLimit)
        );
    }

    #[test]
    fn active_compact_flow_replays_targets_and_final_basis_derivatives() {
        let (context, deformation) = super::super::tests::compact_context();
        let targets = [
            BTreeMap::from([([1, 0], Atom::num(1))]),
            BTreeMap::from([([1, 0], Atom::i())]),
        ];
        let outcome = prepare_weighted_system(
            &context,
            &targets,
            &deformation,
            options(),
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Closed(closed) = outcome else {
            panic!("{outcome:?}")
        };
        closed.differential_system().unwrap().validate().unwrap();
        assert!(closed.diagnostics.rounds > 1);
        for (label, coefficient) in &closed.reduced.targets[0] {
            assert!(
                (&closed.reduced.targets[1][label] - Atom::i() * coefficient)
                    .together()
                    .cancel()
                    .is_zero()
            );
        }
        let decoded = context
            .decode(
                &closed.program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        for target in &targets {
            let result = reduce_sum(
                &decoded,
                target,
                &deformation,
                Default::default(),
                &mut Default::default(),
                &mut BTreeMap::new(),
                &RunContext::default(),
            )
            .unwrap();
            assert!(result.leaves.is_empty() && result.failures.is_empty());
        }
        for basis in &closed.reduced.basis {
            let label = std::array::from_fn(|axis| i64::from(basis.0[axis]));
            let result = reduce_sum(
                &decoded,
                &deformation.derivative(label).unwrap(),
                &deformation,
                Default::default(),
                &mut Default::default(),
                &mut BTreeMap::new(),
                &RunContext::default(),
            )
            .unwrap();
            assert!(result.leaves.is_empty() && result.failures.is_empty());
        }
    }

    #[test]
    fn active_work_caps_leave_an_explicit_unclosed_result() {
        let (context, deformation) = super::super::tests::compact_context();
        for options in [
            WeightedClosureOptions {
                max_rounds: 1,
                ..options()
            },
            WeightedClosureOptions {
                max_requested: 1,
                ..options()
            },
            WeightedClosureOptions {
                application: GuardedReductionLimits {
                    max_rule_applications: 0,
                    ..Default::default()
                },
                ..options()
            },
        ] {
            let result = prepare_weighted_system(
                &context,
                &[BTreeMap::from([([1, 0], Atom::num(1))])],
                &deformation,
                options,
                &RunContext::default(),
            )
            .unwrap();
            assert!(matches!(result, WeightedClosureOutcome::Unresolved(_)));
        }
    }

    #[test]
    fn cancelled_target_retires_uncovered_individual_rows_without_masters() {
        use crate::finite_density::guarded::{
            GuardedIdentity, GuardedIdentityTerm, GuardedMeasureIdentity,
        };
        let eta = symbol!("active_cancel_eta");
        let x = symbol!("active_cancel_x");
        let roles = [IndexRole::Ordinary; 2];
        let context = GuardedContext::new(
            GuardedMeasureIdentity {
                measure: "synthetic pair with equal moments".into(),
                support: "symbolic guarded fixture".into(),
                orientation: "unit".into(),
                normalization: "unit".into(),
                branch: "x nonzero".into(),
                deformation: "fixed difference is zero for every eta".into(),
            },
            roles,
            [symbol!("active_cancel_a"), symbol!("active_cancel_b")],
            vec![eta, x],
            vec![GuardedIdentity {
                id: "equal-moments".into(),
                domain: IndexDomain::new([IndexBounds::fixed(1), IndexBounds::fixed(0)]).unwrap(),
                terms: vec![
                    GuardedIdentityTerm {
                        shift: [0, 0],
                        coefficient: Atom::num(1),
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
            AuxiliaryConvention::EuclideanPlusT,
            roles,
            [true, false],
        )
        .unwrap();
        let directory = std::env::temp_dir().join(format!(
            "rustflow-active-retired-test-{}",
            std::process::id()
        ));
        let outcome = prepare_weighted_system(
            &context,
            &[BTreeMap::from([
                ([1, 0], Atom::num(1)),
                ([0, 1], Atom::num(-1)),
            ])],
            &deformation,
            WeightedClosureOptions {
                checkpoints: Some(directory.clone()),
                ..options()
            },
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Closed(closed) = outcome else {
            panic!("{outcome:?}")
        };
        assert!(closed.reduced.basis.is_empty());
        assert!(closed.reduced.targets[0].is_empty());
        assert!(closed.reduced.candidates.is_empty());
        assert!(closed.differential_system().is_none());
        assert_eq!(closed.diagnostics.retired_unresolved, 2);
        assert!(
            closed
                .reduced
                .nonzero_conditions
                .iter()
                .any(|condition| condition.contains_symbol(x))
        );
        let retired: serde_json::Value = serde_json::from_slice(
            &std::fs::read(directory.join("retired-unresolved.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(retired.as_array().unwrap().len(), 2);
        assert!(retired.as_array().unwrap().iter().all(|entry| {
            entry["unresolved"]
                .as_array()
                .unwrap()
                .iter()
                .all(|r| r["reason"] == "NoApplicableRule")
        }));
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn skipped_cancelled_target_history_retains_required_guard_and_explicit_unassessed_labels() {
        use crate::finite_density::guarded::{
            GuardedIdentity, GuardedIdentityTerm, GuardedMeasureIdentity,
        };
        let eta = symbol!("active_cancel_eta");
        let x = symbol!("active_cancel_x");
        let roles = [IndexRole::Ordinary; 2];
        let context = GuardedContext::new(
            GuardedMeasureIdentity {
                measure: "synthetic pair with equal moments".into(),
                support: "symbolic guarded fixture".into(),
                orientation: "unit".into(),
                normalization: "unit".into(),
                branch: "x nonzero".into(),
                deformation: "fixed difference is zero for every eta".into(),
            },
            roles,
            [symbol!("active_cancel_a"), symbol!("active_cancel_b")],
            vec![eta, x],
            vec![GuardedIdentity {
                id: "equal-moments".into(),
                domain: IndexDomain::new([IndexBounds::fixed(1), IndexBounds::fixed(0)]).unwrap(),
                terms: vec![
                    GuardedIdentityTerm {
                        shift: [0, 0],
                        coefficient: Atom::num(1),
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
            AuxiliaryConvention::EuclideanPlusT,
            roles,
            [true, false],
        )
        .unwrap();
        let directory = std::env::temp_dir().join(format!(
            "rustflow-active-unassessed-test-{}",
            std::process::id()
        ));
        let outcome = prepare_weighted_system(
            &context,
            &[BTreeMap::from([
                ([1, 0], Atom::num(1)),
                ([0, 1], Atom::num(-1)),
            ])],
            &deformation,
            WeightedClosureOptions {
                checkpoints: Some(directory.clone()),
                max_history_candidate_maps: Some(0),
                ..options()
            },
            &RunContext::default(),
        )
        .unwrap();
        let WeightedClosureOutcome::Closed(closed) = outcome else {
            panic!("{outcome:?}")
        };
        assert!(closed.reduced.basis.is_empty());
        assert!(closed.reduced.targets[0].is_empty());
        assert!(closed.reduced.candidates.is_empty());
        assert!(closed.differential_system().is_none());
        assert_eq!(closed.diagnostics.retired_unresolved, 0);
        assert_eq!(closed.diagnostics.history_candidate_assessed, 0);
        assert_eq!(closed.diagnostics.history_candidate_unassessed, 2);
        assert!(
            closed
                .reduced
                .nonzero_conditions
                .iter()
                .any(|condition| condition.contains_symbol(x))
        );
        let retired: serde_json::Value = serde_json::from_slice(
            &std::fs::read(directory.join("retired-unresolved.json")).unwrap(),
        )
        .unwrap();
        assert!(retired.as_array().unwrap().is_empty());
        let metadata: serde_json::Value = serde_json::from_slice(
            &std::fs::read(directory.join(format!(
                "round-{:03}-closed.json",
                closed.diagnostics.rounds - 1
            )))
            .unwrap(),
        )
        .unwrap();
        assert_eq!(metadata["max_history_candidate_maps"], 0);
        let coverage = &metadata["active_state"]["historical_candidate_collection"];
        assert!(coverage["selected"].as_array().unwrap().is_empty());
        assert_eq!(coverage["unassessed"], serde_json::json!([[0, 1], [1, 0]]));
        assert_eq!(coverage["unassessed_are_unresolved"], false);
        assert_eq!(coverage["history_is_coverage_certificate"], false);
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn direct_zero_prepass_supersedes_a_retained_recurrence_and_memoizes_search() {
        use crate::finite_density::guarded::{
            GuardedIdentity, GuardedIdentityTerm, GuardedMeasureIdentity,
        };
        let eta = symbol!("active_zero_eta");
        let a = symbol!("active_zero_a");
        let roles = [IndexRole::Ordinary];
        let context = GuardedContext::new(
            GuardedMeasureIdentity {
                measure: "synthetic zero sequence".into(),
                support: "positive indices".into(),
                orientation: "unit".into(),
                normalization: "unit".into(),
                branch: "eta nonzero".into(),
                deformation: "single shifted ordinary factor".into(),
            },
            roles,
            [a],
            vec![eta],
            vec![
                GuardedIdentity {
                    id: "early-recurrence".into(),
                    domain: IndexDomain::new([IndexBounds::new(Some(1), None).unwrap()]).unwrap(),
                    terms: vec![
                        GuardedIdentityTerm {
                            shift: [0],
                            coefficient: Atom::num(1),
                        },
                        GuardedIdentityTerm {
                            shift: [-1],
                            coefficient: Atom::num(-1),
                        },
                    ],
                    nonzero_conditions: vec![],
                },
                GuardedIdentity {
                    id: "later-zero".into(),
                    domain: IndexDomain::new([IndexBounds::fixed(1)]).unwrap(),
                    terms: vec![GuardedIdentityTerm {
                        shift: [1],
                        coefficient: Atom::var(eta),
                    }],
                    nonzero_conditions: vec![],
                },
            ],
        )
        .unwrap();
        let deformation =
            FixedShellDeformation::new(eta, AuxiliaryConvention::EuclideanPlusT, roles, [true])
                .unwrap();
        let domains = [IndexDomain::new([IndexBounds::fixed(2)]).unwrap()];
        let old = context
            .discover(domains.to_vec(), [], Default::default())
            .unwrap()
            .program;
        assert!(
            !old.reduce([2], Default::default())
                .unwrap()
                .unresolved
                .is_empty()
        );
        let options = WeightedClosureOptions {
            max_direct_zero_attempts: 64,
            ..options()
        };
        let requested = BTreeSet::from([[2]]);
        let mut retained = Some(old);
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let mut conditions = BTreeMap::new();
        let mut learned = Vec::new();
        let mut gaps = BTreeMap::new();
        let found = discover_with_refinement(
            &context,
            &domains,
            &[],
            &requested,
            &deformation,
            &options,
            &mut learned,
            &mut retained,
            &mut gaps,
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        let reduced = found.program.reduce([2], Default::default()).unwrap();
        assert!(reduced.terms.is_empty() && reduced.unresolved.is_empty());
        assert!(
            reduced
                .nonzero_conditions
                .iter()
                .any(|condition| condition.contains_symbol(eta))
        );
        assert_eq!(diagnostics.direct_zero_rules, 1);
        assert!(diagnostics.direct_zero_completed_points.contains(&vec![2]));
        let attempts = diagnostics.direct_zero_attempts;
        retained = Some(found.program);
        let repeated = discover_with_refinement(
            &context,
            &domains,
            &[],
            &requested,
            &deformation,
            &options,
            &mut learned,
            &mut retained,
            &mut gaps,
            &mut diagnostics,
            &mut conditions,
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(diagnostics.direct_zero_attempts, attempts);
        assert!(
            repeated
                .program
                .reduce([2], Default::default())
                .unwrap()
                .unresolved
                .is_empty()
        );
    }
}
