//! Bounded native ray discovery with an exact-point fallback for each request.
//!
//! The retained program is frozen throughout a batch. Every fresh proof is
//! replayed against the original source context; no row elimination lives here.
use super::*;
use rustred::solver::guarded::GuardedApplicationStatus;
use serde_json::{Value, json};

fn add(counter: &mut usize, amount: usize) -> Result<()> {
    *counter = counter
        .checked_add(amount)
        .ok_or_else(|| Error::Limit("requested discovery counter overflow".into()))?;
    Ok(())
}

fn record_gaps<const N: usize>(
    historical: &mut HistoricalDiscoveryGaps<N>,
    gaps: Vec<GuardedUnresolved<N>>,
) -> Vec<Value> {
    gaps.into_iter()
        .map(|gap| {
            let bounds = std::array::from_fn(|axis| {
                let bound = gap.domain.bounds()[axis];
                [bound.lower(), bound.upper()]
            });
            let record = json!({"bounds":bounds.to_vec(),
                "reason":format!("{:?}",gap.reason), "detail":gap.detail});
            historical.insert((bounds, format!("{:?}:{}", gap.reason, gap.detail)), gap);
            record
        })
        .collect()
}

/// Return true only for the two failures for which this schedule permits a
/// fresh original-source search. A successful application is not closure.
fn probe<const N: usize>(
    program: &GuardedReductionProgram<N>,
    point: [i64; N],
    deformation: &FixedShellDeformation<N>,
    conditions: &mut BTreeMap<String, Atom>,
    run: &RunContext,
    record: &mut Value,
) -> Result<bool> {
    run.cancellation.check()?;
    let applied = program.apply_for_discovery(point)?;
    *record = json!({"status":format!("{:?}",applied.status),
        "rhs":applied.terms.iter().map(|(i,c)|json!({"integral":i.to_vec(),
            "coefficient":c.to_canonical_string()})).collect::<Vec<_>>(),
        "nonzero_conditions":applied.nonzero_conditions.iter()
            .map(|condition| condition.to_canonical_string()).collect::<Vec<_>>()});
    for label in applied.terms.keys() {
        deformation.validate_integral(label)?;
    }
    // Failed native applications return no asserted conditions. Keep every
    // condition from successful probes, including ones whose RHS later cancels.
    retain_conditions(conditions, applied.nonzero_conditions);
    match applied.status {
        GuardedApplicationStatus::Applied { .. } | GuardedApplicationStatus::Zero => Ok(false),
        GuardedApplicationStatus::Unresolved(
            GuardedApplicationFailure::NoApplicableRule
            | GuardedApplicationFailure::ConditionVanished { .. },
        ) => Ok(true),
        status => Err(Error::IncompleteReduction(format!(
            "native requested-point application failed at {point:?}: {status:?}"
        ))),
    }
}

/// The caller runs the shared direct-zero prepass first. This helper keeps its
/// zero-first/retained-first order and returns completed versus deferred point
/// transactions separately. Neither set is a coverage or master certificate.
#[allow(clippy::too_many_arguments)]
pub(super) fn prepare<const N: usize>(
    context: &GuardedContext<N>,
    requested: &BTreeSet<[i64; N]>,
    deformation: &FixedShellDeformation<N>,
    options: &WeightedClosureOptions,
    retained: &mut Option<GuardedReductionProgram<N>>,
    gaps: &mut HistoricalDiscoveryGaps<N>,
    diagnostics: &mut WeightedClosureDiagnostics,
    conditions: &mut BTreeMap<String, Atom>,
    run: &RunContext,
) -> Result<(GuardedDiscovery<N>, BTreeSet<[i64; N]>, BTreeSet<[i64; N]>)> {
    let policy = options.requested_ray_point.ok_or_else(|| {
        Error::InvalidInput("requested ray/point discovery requires its explicit options".into())
    })?;
    // Admission is independent of budget and lexical position. A zero budget
    // must not hide an invalid later request or a nonzero storage tail.
    for point in requested {
        deformation.validate_integral(point)?;
    }
    let reservation = policy
        .max_domains_per_ray
        .checked_add(policy.max_domains_per_point)
        .ok_or_else(|| Error::Limit("requested ray/point reservation overflow".into()))?;
    if policy.max_domains_per_ray == 0 || policy.max_domains_per_point == 0 {
        return Err(Error::InvalidInput(
            "requested ray/point discovery needs positive per-search caps".into(),
        ));
    }
    add(&mut diagnostics.requested_ray_point_calls, 1)?;
    let call = diagnostics.requested_ray_point_calls;
    let mut remaining = options.discovery.max_domains;
    let mut transactions = Vec::new();
    let mut completed = BTreeSet::new();
    let mut deferred = BTreeSet::new();
    let mut result = (|| -> Result<GuardedDiscovery<N>> {
        run.cancellation.check()?;
        let frozen = match retained.take() {
            Some(program)
                if std::sync::Arc::ptr_eq(program.native().sources(), context.sources())
                    && program.native().terminals().is_empty() =>
            {
                program
            }
            Some(program) => {
                // Bind to the caller's complete source context even if every
                // request is already covered and no fresh rules are found.
                // Empty right-hand rules preserve the retained order, and the
                // explicit empty stopping set removes stale terminals.
                run.cancellation.check()?;
                add(&mut diagnostics.native_rule_unions, 1)?;
                let empty = context.discover(Vec::new(), [], options.discovery)?.program;
                program.union_verified(empty, [], options.max_reused_rules)?
            }
            None => context.discover(Vec::new(), [], options.discovery)?.program,
        };
        let mut total_rules = frozen.native().rules().len();
        if total_rules > options.max_reused_rules {
            return Err(Error::Limit(
                "requested native program exceeds the retained-rule budget".into(),
            ));
        }
        let mut fresh = Vec::new();
        for &point in requested {
            let mut transaction = json!({"integral":point.to_vec(), "complete":false,
                "allocated_ray_domains":0, "allocated_point_domains":0,
                "remaining_before":remaining});
            let attempted = (|| -> Result<()> {
                if !probe(
                    &frozen,
                    point,
                    deformation,
                    conditions,
                    run,
                    &mut transaction["retained_probe"],
                )? {
                    transaction["complete"] = json!(true);
                    transaction["reason"] = json!("retained rule or declared zero applies");
                    completed.insert(point);
                    return Ok(());
                }
                if remaining < reservation {
                    transaction["reason"] = json!("ray/point reservation exceeds remaining budget");
                    deferred.insert(point);
                    add(&mut diagnostics.requested_ray_point_deferred, 1)?;
                    return Ok(());
                }
                let ray =
                    discovery_domains(&BTreeSet::from([point]), &deformation.roles, true, true)?
                        .pop()
                        .expect("one admitted request produces one ray")
                        .intersection(deformation.admitted_domain())
                        .ok_or_else(|| {
                            Error::InvalidInput("requested ray misses its admitted point".into())
                        })?;
                run.cancellation.check()?;
                remaining -= policy.max_domains_per_ray;
                add(
                    &mut diagnostics.requested_ray_point_allocated_domains,
                    policy.max_domains_per_ray,
                )?;
                transaction["allocated_ray_domains"] = json!(policy.max_domains_per_ray);
                let ray_found = context.discover(
                    vec![ray],
                    [],
                    GuardedDiscoveryOptions {
                        max_domains: policy.max_domains_per_ray,
                        ..options.discovery
                    },
                )?;
                transaction["ray_gaps"] = json!(record_gaps(gaps, ray_found.unresolved));
                let fallback = probe(
                    &ray_found.program,
                    point,
                    deformation,
                    conditions,
                    run,
                    &mut transaction["ray_probe"],
                )?;
                add(&mut total_rules, ray_found.program.native().rules().len())?;
                if total_rules > options.max_reused_rules {
                    return Err(Error::Limit(
                        "requested native program exceeds the retained-rule budget".into(),
                    ));
                }
                fresh.push(ray_found.program);
                if fallback {
                    run.cancellation.check()?;
                    remaining -= policy.max_domains_per_point;
                    add(
                        &mut diagnostics.requested_ray_point_allocated_domains,
                        policy.max_domains_per_point,
                    )?;
                    add(&mut diagnostics.requested_ray_point_fallbacks, 1)?;
                    transaction["allocated_point_domains"] = json!(policy.max_domains_per_point);
                    let point_found = context.discover(
                        vec![
                            IndexDomain::new(point.map(IndexBounds::fixed))
                                .map_err(|error| Error::InvalidInput(error.to_string()))?,
                        ],
                        [],
                        GuardedDiscoveryOptions {
                            max_domains: policy.max_domains_per_point,
                            ..options.discovery
                        },
                    )?;
                    transaction["point_gaps"] = json!(record_gaps(gaps, point_found.unresolved));
                    let still_unresolved = probe(
                        &point_found.program,
                        point,
                        deformation,
                        conditions,
                        run,
                        &mut transaction["point_probe"],
                    )?;
                    transaction["point_still_unresolved"] = json!(still_unresolved);
                    add(&mut total_rules, point_found.program.native().rules().len())?;
                    if total_rules > options.max_reused_rules {
                        return Err(Error::Limit(
                            "requested native program exceeds the retained-rule budget".into(),
                        ));
                    }
                    fresh.push(point_found.program);
                }
                completed.insert(point);
                transaction["complete"] = json!(true);
                transaction["reason"] =
                    json!("bounded native search completed; no coverage inferred");
                Ok(())
            })();
            transaction["remaining_after"] = json!(remaining);
            if let Err(error) = &attempted {
                transaction["error"] = json!(error.to_string());
            }
            transactions.push(transaction);
            attempted?;
        }
        // Pairwise union preserves ray-before-point and request order while
        // keeping vector concatenation balanced. Leaf programs have already
        // passed native source replay; verified composition retains that proof.
        while fresh.len() > 1 {
            let mut next = Vec::new();
            let mut iter = fresh.into_iter();
            while let Some(left) = iter.next() {
                next.push(if let Some(right) = iter.next() {
                    run.cancellation.check()?;
                    add(&mut diagnostics.native_rule_unions, 1)?;
                    left.union_verified(right, [], options.max_reused_rules)?
                } else {
                    left
                });
            }
            fresh = next;
        }
        run.cancellation.check()?;
        let program = if let Some(fresh) = fresh.pop() {
            add(&mut diagnostics.native_rule_unions, 1)?;
            frozen.union_verified(fresh, [], options.max_reused_rules)?
        } else {
            frozen
        };
        if program.native().rules().len() > options.max_reused_rules {
            return Err(Error::Limit(
                "requested native program exceeds the retained-rule budget".into(),
            ));
        }
        run.cancellation.check()?;
        Ok(GuardedDiscovery {
            program,
            unresolved: gaps.values().cloned().collect(),
        })
    })();
    diagnostics.historical_discovery_domains = gaps.len();
    // An aborted call gets an explicitly incomplete evidence record, never a
    // completed-round marker. The active owner publishes the program-bound
    // round checkpoint after receiving this complete result.
    if let Some(directory) = &options.checkpoints {
        std::fs::create_dir_all(directory)?;
        let mut program_binding = result
            .as_ref()
            .ok()
            .map(|found| {
                let bytes = found.program.encode(Default::default())?;
                Ok::<_, Error>(json!({"blake3":blake3::hash(&bytes).to_hex().to_string(),
                "bytes":bytes.len(), "rules":found.program.native().rules().len()}))
            })
            .transpose()?;
        if result.is_ok() {
            if let Err(error) = run.cancellation.check() {
                result = Err(error);
                program_binding = None;
            }
        }
        let metadata = json!({"schema":1, "schedule":"requested-ray-point-v1", "call":call,
            "source_measure_id":context.sources().measure_id(),
            "physical_arity":context.physical_arity(), "storage_capacity":N,
            "policy":policy, "discovery":options.discovery,
            "max_reused_rules":options.max_reused_rules,
            "status":if result.is_ok() { "complete" } else { "incomplete" },
            "error":result.as_ref().err().map(ToString::to_string),
            "native_program":program_binding, "empty_terminals":true,
            "rule_precedence":"shared direct zeros, frozen retained rules, fresh request order (ray then point)",
            "allocated_domains":options.discovery.max_domains-remaining,
            "remaining_domains":remaining, "allocation_is_measured_visits":false,
            "completed_requests":completed.iter().map(|i|i.to_vec()).collect::<Vec<_>>(),
            "deferred_requests":deferred.iter().map(|i|i.to_vec()).collect::<Vec<_>>(),
            "transactions":transactions, "historical_gap_count":gaps.len(),
            "completed_search_implies_coverage":false});
        let path = directory.join(format!("requested-discovery-{call:04}.json"));
        let temporary = path.with_extension("json.part");
        std::fs::write(
            &temporary,
            serde_json::to_vec_pretty(&metadata)
                .map_err(|error| Error::Cache(error.to_string()))?,
        )?;
        std::fs::rename(temporary, path)?;
    }
    result.map(|found| (found, completed, deferred))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_density::guarded::{
        GuardedIdentity, GuardedIdentityTerm, GuardedMeasureIdentity,
    };

    #[test]
    fn original_two_row_corpus_needs_exact_point_after_bounded_ray_miss() {
        // This deletion-minimal two-row excerpt was recovered from the real
        // singleton pilot. The Ward and Lorentz rows together give I=I_child
        // at this point; neither source alone yields a bounded point rule.
        // All original guards and declared zero boxes remain unchanged.
        let fixture: Value = serde_json::from_str(include_str!(
            "../../../fixtures/finite_density/requested_point_source_rows.json"
        ))
        .unwrap();
        let domain = |value: &Value| {
            let bounds: Vec<[Option<i64>; 2]> = serde_json::from_value(value.clone()).unwrap();
            IndexDomain::<12>::new(std::array::from_fn(|i| {
                IndexBounds::new(bounds[i][0], bounds[i][1]).unwrap()
            }))
            .unwrap()
        };
        let parse =
            |text: &str| Atom::parse(text, "requested_point_fixture", Default::default()).unwrap();
        let epsilon = symbol!("requested_point_fixture::epsilon");
        let eta = symbol!("requested_point_fixture::eta");
        let indices =
            std::array::from_fn(|axis| symbol!(format!("requested_point_fixture::a_{axis}")));
        let roles = std::array::from_fn(|axis| match fixture["roles"][axis].as_u64().unwrap() {
            0 => IndexRole::Ordinary,
            1 => IndexRole::RequiredCut,
            2 => IndexRole::Occupation,
            _ => panic!("invalid fixture role"),
        });
        let sources = fixture["source_rows"]
            .as_array()
            .unwrap()
            .iter()
            .map(|row| GuardedIdentity {
                id: row["source_id"].as_str().unwrap().into(),
                domain: domain(&row["domain"]),
                nonzero_conditions: row["conditions"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|condition| parse(condition.as_str().unwrap()))
                    .collect(),
                terms: row["terms"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|term| GuardedIdentityTerm {
                        shift: std::array::from_fn(|axis| {
                            assert_eq!(term["label"][axis][0], true);
                            term["label"][axis][1].as_i64().unwrap().try_into().unwrap()
                        }),
                        coefficient: parse(term["coefficient"].as_str().unwrap()),
                    })
                    .collect(),
            })
            .collect();
        let context = GuardedContext::new_with_physical_arity(
            GuardedMeasureIdentity {
                measure: "original two-row polynomial singleton source excerpt".into(),
                support: "original lower-contact and free-virtual zero boxes".into(),
                orientation: "original future massless singleton".into(),
                normalization: "homogeneous rows; no period values".into(),
                branch: "original sealed high-dimensional source continuation".into(),
                deformation: "original four uncut physical factors D-eta".into(),
            },
            roles,
            indices,
            vec![epsilon, eta],
            sources,
            11,
        )
        .unwrap()
        .with_measure_zero_domains(
            fixture["zero_domains"]
                .as_array()
                .unwrap()
                .iter()
                .map(domain)
                .collect(),
        )
        .unwrap();
        let mut bounds = *IndexDomain::for_roles(&roles).bounds();
        for bound in &mut bounds[5..9] {
            *bound = IndexBounds::new(None, Some(0)).unwrap();
        }
        bounds[11] = IndexBounds::fixed(0);
        let deformation = FixedShellDeformation::new(
            eta,
            AuxiliaryConvention::NativeMinusEta,
            roles,
            std::array::from_fn(|axis| (1..=4).contains(&axis)),
        )
        .unwrap()
        .with_physical_arity(11)
        .unwrap()
        .with_admitted_domain(IndexDomain::new(bounds).unwrap())
        .unwrap();
        let point = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0];
        let requested = BTreeSet::from([point]);
        let ray = discovery_domains(&requested, &roles, true, true).unwrap();
        let ray_only = context
            .discover(
                ray,
                [],
                GuardedDiscoveryOptions {
                    max_depth: 3,
                    max_domains: 1,
                    sample_seed: 0,
                },
            )
            .unwrap();
        assert!(ray_only.program.native().rules().is_empty());
        assert!(matches!(
            ray_only.program.apply_for_discovery(point).unwrap().status,
            GuardedApplicationStatus::Unresolved(GuardedApplicationFailure::NoApplicableRule)
        ));
        assert!(!ray_only.unresolved.is_empty());
        let mut diagnostics = WeightedClosureDiagnostics::default();
        let (found, completed, deferred) = prepare(
            &context,
            &requested,
            &deformation,
            &WeightedClosureOptions {
                requested_ray_point: Some(RequestedRayPointOptions {
                    max_domains_per_ray: 1,
                    max_domains_per_point: 1,
                }),
                discovery: GuardedDiscoveryOptions {
                    max_depth: 3,
                    max_domains: 2,
                    sample_seed: 0,
                },
                max_reused_rules: 16,
                ..Default::default()
            },
            &mut None,
            &mut BTreeMap::new(),
            &mut diagnostics,
            &mut BTreeMap::new(),
            &RunContext::default(),
        )
        .unwrap();
        assert_eq!(completed, requested);
        assert!(deferred.is_empty());
        assert_eq!(diagnostics.requested_ray_point_allocated_domains, 2);
        assert_eq!(diagnostics.requested_ray_point_fallbacks, 1);
        assert_eq!(found.program.native().rules().len(), 1);
        assert_eq!(found.unresolved.len(), ray_only.unresolved.len());
        let decoded = context
            .decode(
                &found.program.encode(Default::default()).unwrap(),
                Default::default(),
            )
            .unwrap();
        let application = decoded.apply_for_discovery(point).unwrap();
        assert!(matches!(
            application.status,
            GuardedApplicationStatus::Applied { .. }
        ));
        assert_eq!(
            application.terms,
            BTreeMap::from([([1, 1, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0], Atom::one())])
        );
        assert!(application.nonzero_conditions.is_empty());
        assert!(decoded.native().terminals().is_empty());
    }
}
