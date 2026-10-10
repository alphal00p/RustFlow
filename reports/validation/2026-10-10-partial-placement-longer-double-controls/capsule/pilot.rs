//! Fresh partial-placement native active-closure diagnostic, no physical flow admission.
use rustred::persistence::BinaryIoLimits;
use rustred::solver::SearchOptions;
use rustred::solver::guarded::{
    GuardedApplicationFailure as Failure, GuardedApplicationStatus as Status, GuardedProgram,
    IndexBounds, IndexDomain,
};
mod prepare;
fn bounds<const N: usize>(v: &IndexDomain<N>) -> Vec<(Option<i64>, Option<i64>)> {
    v.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}
fn fallback(s: &Status) -> bool {
    matches!(
        s,
        Status::Unresolved(Failure::NoApplicableRule | Failure::ConditionVanished { .. })
    )
}
use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

type Point = [i64; 16];
type Row = BTreeMap<Point, Atom>;
struct Reduced {
    leaves: Row,
    terms: Row,
    conditions: BTreeSet<String>,
    failures: Vec<serde_json::Value>,
    applications: usize,
}
fn add(row: &mut Row, i: Point, c: Atom) {
    *row.entry(i).or_default() += c;
}
fn clean(row: &mut Row) {
    row.retain(|_, c| {
        *c = c.together().cancel();
        !c.is_zero()
    });
}
fn reduce(program: &GuardedProgram<16>, row: &Row, admitted: &IndexDomain<16>) -> Reduced {
    let mut r = Reduced {
        leaves: Row::new(),
        terms: Row::new(),
        conditions: BTreeSet::new(),
        failures: Vec::new(),
        applications: 0,
    };
    for (i, c) in row {
        if c.is_zero() {
            continue;
        }
        assert!(admitted.contains(i));
        let x = program.reduce(*i, Default::default()).unwrap();
        r.applications += x.rule_applications;
        r.conditions.extend(
            x.nonzero_conditions
                .iter()
                .map(|c| c.to_expression().to_string()),
        );
        for (j, w) in x.terms {
            assert!(admitted.contains(&j));
            add(&mut r.terms, j, c * w.to_expression());
        }
        for t in x.unresolved {
            assert!(admitted.contains(&t.integral));
            if t.reason == rustred::solver::guarded::GuardedApplicationFailure::NoApplicableRule {
                add(&mut r.leaves, t.integral, c * t.coefficient.to_expression());
            } else {
                r.failures.push(serde_json::json!({"indices":t.integral,"coefficient":(c*t.coefficient.to_expression()).to_string(),"reason":format!("{:?}",t.reason)}));
            }
        }
    }
    clean(&mut r.leaves);
    clean(&mut r.terms);
    r
}
fn derivative(i: Point, shifted: &[usize]) -> Row {
    let mut row = Row::new();
    for &slot in shifted {
        if i[slot] != 0 {
            let mut j = i;
            j[slot] += 1;
            add(&mut row, j, Atom::num(i[slot]));
        }
    }
    clean(&mut row);
    row
}
fn display(row: &Row) -> Vec<serde_json::Value> {
    row.iter()
        .map(|(i, c)| serde_json::json!({"indices":i,"coefficient":c.to_string()}))
        .collect()
}
fn setting(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .map(|s| s.parse().unwrap())
        .unwrap_or(default)
}
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let out = std::path::PathBuf::from(&args[2]);
    assert!(!out.exists());
    std::fs::create_dir_all(&out).unwrap();
    let limits = BinaryIoLimits::default();
    let prepared = prepare::prepare(&args, &out);
    let sources = prepared.sources;
    let targets = prepared.targets;
    let shifted = prepared.shifted;
    let admitted = prepared.admitted;
    let input_path = &args[1];
    let old_proof_replay: Option<usize> = None;
    let max_requested = setting("PILOT_REQUESTS", 4096);
    let max_rounds = setting("PILOT_ROUNDS", 16);
    let max_frontier = setting("PILOT_FRONTIER", 1024);
    let max_rules = setting("PILOT_RULES", 65536);
    let zero_attempts = setting("PILOT_ZERO_ATTEMPTS", 8192);
    let ray_budget = setting("PILOT_RAY_DOMAINS", 1);
    let point_budget = setting("PILOT_POINT_DOMAINS", 1);
    let depth = setting("PILOT_DEPTH", 3) as u32;
    let mut program = GuardedProgram::new(sources.clone(), vec![], []).unwrap();
    let mut frontier = BTreeSet::new();
    let mut all_conditions = BTreeSet::new();
    for row in &targets {
        let red = reduce(&program, row, &admitted);
        assert!(
            red.failures.is_empty(),
            "native failures: {:?}",
            red.failures
        );
        frontier.extend(red.leaves.keys().copied());
    }
    std::fs::write(out.join("configuration.json"),serde_json::to_vec_pretty(&serde_json::json!({"request_history_policy":"explicit requested labels only; basis and raw derivative labels deferred until requested; max_requested caps cumulative history","physical_input":args[1],"source_count":sources.sources().len(),"original_targets":targets.iter().map(display).collect::<Vec<_>>(),"old_proofs_reused":false,"old_proof_optional_replay":old_proof_replay,"max_requested":max_requested,"max_rounds":max_rounds,"max_frontier":max_frontier,"max_rules":max_rules,"zero_attempts_per_round":zero_attempts,"depth":depth,"ray_domains_per_requested_point":ray_budget,"point_domains_per_fallback":point_budget,"shifted_ordinary_slots":shifted,"measure_id":sources.measure_id()})).unwrap()).unwrap();
    let mut history = Vec::new();
    let mut closed = None;
    let mut stop = "round-budget";
    // Exact request history is scheduling evidence, not a coverage certificate.
    // Start with the same three original labels and their twelve raw derivative
    // labels as the preserved first pair; never mark a reduction output searched.
    let roots = targets
        .iter()
        .flat_map(|r| r.keys().copied())
        .collect::<BTreeSet<_>>();
    let mut requested = roots.clone();
    for &i in &roots {
        requested.extend(derivative(i, &shifted).keys().copied());
    }
    let mut attempted_points = BTreeSet::new();
    for round in 0..max_rounds {
        let started = std::time::Instant::now();
        let before = frontier.clone();
        if attempted_points.union(&requested).count() > max_requested {
            stop = "request-budget";
            break;
        }
        let points = requested.iter().copied().collect::<Vec<_>>();
        attempted_points.extend(points.iter().copied());
        let zero = sources
            .direct_zero_rules_at_points(&points, zero_attempts)
            .unwrap();
        let zero_gaps=zero.solution.unresolved.iter().map(|g|serde_json::json!({"domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
        let zero_count = zero.solution.rules.len();
        let zero_ids=zero.solution.rules.iter().map(|r|serde_json::json!({"target":r.candidate().target.powers().iter().map(|p|p.value()).collect::<Vec<_>>(),"sources":r.candidate().sources.iter().map(|s|serde_json::json!({"source_id":sources.sources()[s.basis_row].id,"seed":s.seed.integral.powers().iter().map(|p|p.value()).collect::<Vec<_>>()})).collect::<Vec<_>>()})).collect::<Vec<_>>();
        let zero_program = GuardedProgram::new(sources.clone(), zero.solution.rules, []).unwrap();
        let mut combined = zero_program.union_verified(program, [], max_rules).unwrap();
        let mut fresh_programs = Vec::new();
        let mut probes = Vec::new();
        for point in points {
            let existing = combined.apply(&point).unwrap();
            if !fallback(&existing.status) {
                assert!(matches!(
                    existing.status,
                    Status::Applied { .. } | Status::Zero
                ));
                continue;
            }
            let ray = IndexDomain::new(std::array::from_fn(|i| {
                if point[i] == 0 {
                    IndexBounds::fixed(0)
                } else if point[i] > 0 {
                    IndexBounds::new(Some(point[i]), None).unwrap()
                } else {
                    IndexBounds::new(None, Some(point[i])).unwrap()
                }
            }))
            .unwrap()
            .intersection(&admitted)
            .unwrap();
            let found = sources
                .solve_domains(
                    vec![ray],
                    SearchOptions {
                        max_depth: Some(depth),
                        sample_seed: 0,
                        ..Default::default()
                    },
                    ray_budget,
                )
                .unwrap();
            let mut gaps=found.unresolved.iter().map(|g|serde_json::json!({"phase":"ray","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})).collect::<Vec<_>>();
            let ray_program = GuardedProgram::new(sources.clone(), found.rules, []).unwrap();
            let ray_status = ray_program.apply(&point).unwrap().status;
            let fallback = fallback(&ray_status);
            assert!(fallback || matches!(ray_status, Status::Applied { .. } | Status::Zero));
            fresh_programs.push(ray_program);
            let mut point_status = None;
            let portfolio_stats: Vec<serde_json::Value> = Vec::new();
            if fallback {
                let found = sources
                    .solve_domains(
                        vec![IndexDomain::new(point.map(IndexBounds::fixed)).unwrap()],
                        SearchOptions {
                            max_depth: Some(depth),
                            sample_seed: 0,
                            ..Default::default()
                        },
                        point_budget,
                    )
                    .unwrap();
                gaps.extend(found.unresolved.iter().map(|g|serde_json::json!({"phase":"point","domain":bounds(&g.domain),"reason":format!("{:?}",g.reason),"detail":g.detail})));
                let exact = GuardedProgram::new(sources.clone(), found.rules, []).unwrap();
                point_status = Some(format!("{:?}", exact.apply(&point).unwrap().status));
                fresh_programs.push(exact);
            }
            probes.push(serde_json::json!({"target":point,"ray_status":format!("{:?}",ray_status),"point_fallback":fallback,"point_status":point_status,"portfolio":portfolio_stats,"gaps":gaps}));
        }
        let fresh_count = fresh_programs
            .iter()
            .map(|p| p.rules().len())
            .sum::<usize>();
        // GuardedRule intentionally has no Clone/export escape hatch. Combine the
        // independently replayed point programs in a balanced native union tree,
        // retaining source/order checks without quadratic sequential replay cost.
        while fresh_programs.len() > 1 {
            let mut next = Vec::new();
            let mut iter = fresh_programs.into_iter();
            while let Some(left) = iter.next() {
                next.push(if let Some(right) = iter.next() {
                    left.union_verified(right, [], max_rules).unwrap()
                } else {
                    left
                });
            }
            fresh_programs = next;
        }
        if let Some(fresh) = fresh_programs.pop() {
            combined = combined.union_verified(fresh, [], max_rules).unwrap();
        }
        let mut retired = Vec::new();
        for &i in &before {
            let red = reduce(&combined, &BTreeMap::from([(i, Atom::num(1))]), &admitted);
            assert!(
                red.failures.is_empty(),
                "native failures: {:?}",
                red.failures
            );
            if !red.leaves.contains_key(&i) {
                retired.push(i);
            }
            all_conditions.extend(red.conditions);
        }
        let mut reachable = BTreeSet::new();
        let mut target_rows = Vec::new();
        for row in &targets {
            let red = reduce(&combined, row, &admitted);
            assert!(
                red.failures.is_empty(),
                "native failures: {:?}",
                red.failures
            );
            reachable.extend(red.leaves.keys().copied());
            all_conditions.extend(red.conditions);
            target_rows.push(display(&red.leaves));
        }
        // Rebuild only original-reachable obligations. Each prospective basis label
        // and each raw derivative label must first have been an explicit native
        // request. Containment in a searched ray is not a substitute for that check.
        frontier.clear();
        let mut needed = roots.clone();
        let mut processed = BTreeSet::new();
        let mut pending_new = BTreeSet::new();
        let mut deferred_basis = BTreeSet::new();
        let mut deferred_raw_derivatives = BTreeSet::new();
        let mut derivative_rows = Vec::new();
        let mut overflow = false;
        while let Some(i) = reachable.pop_first() {
            if !frontier.insert(i) {
                continue;
            }
            needed.insert(i);
            if frontier.len() > max_frontier {
                overflow = true;
                break;
            }
            if !attempted_points.contains(&i) {
                pending_new.insert(i);
                deferred_basis.insert(i);
                continue;
            }
            let raw = derivative(i, &shifted);
            needed.extend(raw.keys().copied());
            let missing = raw
                .keys()
                .filter(|j| !attempted_points.contains(*j))
                .copied()
                .collect::<Vec<_>>();
            if !missing.is_empty() {
                pending_new.extend(missing.iter().copied());
                deferred_raw_derivatives.extend(missing);
                continue;
            }
            processed.insert(i);
            let red = reduce(&combined, &raw, &admitted);
            assert!(
                red.failures.is_empty(),
                "native failures: {:?}",
                red.failures
            );
            all_conditions.extend(red.conditions);
            reachable.extend(
                red.leaves
                    .keys()
                    .filter(|j| !frontier.contains(*j))
                    .copied(),
            );
            derivative_rows.push(serde_json::json!({"basis":i,"raw_derivative":display(&raw),"reduced_actual_derivative":display(&red.leaves)}));
        }
        let report = serde_json::json!({"round":round,"seconds":started.elapsed().as_secs_f64(),"requested":requested,"attempted_exact_points":attempted_points,"next_needed":needed,"deferred_basis_labels":deferred_basis,"deferred_raw_derivative_labels":deferred_raw_derivatives,"source_zero_gaps":zero_gaps,"source_zero_attempts":zero.attempted_rows,"zero_completed_points":zero.completed_points,"zero_rules":zero_ids,"zero_first_rule_count":zero_count,"fresh_rules":fresh_count,"rule_count":combined.rules().len(),"retired_reducible_frontier":retired,"frontier":frontier,"processed_derivatives":processed,"pending_new_derivative_labels":pending_new,"original_target_rows":target_rows,"derivative_rows":derivative_rows,"discovery":probes,"conditions_encountered":all_conditions,"frontier_budget_exhausted":overflow});
        std::fs::write(
            out.join(format!("round-{round:03}.json")),
            serde_json::to_vec_pretty(&report).unwrap(),
        )
        .unwrap();
        std::fs::write(
            out.join(format!("round-{round:03}.bin")),
            combined.encode_native(limits).unwrap(),
        )
        .unwrap();
        history.push(serde_json::json!({"round":round,"frontier":frontier.len(),"rules":combined.rules().len(),"requested":requested.len(),"zero_rules":zero_count,"point_fallbacks":report["discovery"].as_array().unwrap().iter().filter(|p|p["point_fallback"]==true).count(),"retired":report["retired_reducible_frontier"].as_array().unwrap().len(),"seconds":started.elapsed().as_secs_f64()}));
        std::fs::write(
            out.join("progress.json"),
            serde_json::to_vec_pretty(&history).unwrap(),
        )
        .unwrap();
        eprintln!(
            "round={round} frontier={} rules={} requested={} pending={}",
            frontier.len(),
            combined.rules().len(),
            requested.len(),
            pending_new.len()
        );
        if overflow {
            stop = "frontier-budget";
            break;
        }
        if pending_new.is_empty() {
            let audited = combined
                .with_terminals_replayed(frontier.iter().copied(), max_rules)
                .unwrap();
            let mut audits = Vec::new();
            for row in targets
                .iter()
                .cloned()
                .chain(frontier.iter().map(|&i| derivative(i, &shifted)))
            {
                let red = reduce(&audited, &row, &admitted);
                assert!(red.failures.is_empty() && red.leaves.is_empty());
                all_conditions.extend(red.conditions);
                audits.push(display(&red.terms));
            }
            let bytes = audited.encode_native(limits).unwrap();
            std::fs::write(out.join("closed.bin"), &bytes).unwrap();
            let decoded =
                GuardedProgram::decode_generated(&bytes, sources.clone(), limits).unwrap();
            for row in targets
                .iter()
                .cloned()
                .chain(frontier.iter().map(|&i| derivative(i, &shifted)))
            {
                let a = reduce(&audited, &row, &admitted);
                let b = reduce(&decoded, &row, &admitted);
                assert!(
                    a.failures.is_empty()
                        && b.failures.is_empty()
                        && a.leaves.is_empty()
                        && b.leaves.is_empty()
                );
                assert_eq!(a.terms, b.terms);
                assert_eq!(a.conditions, b.conditions);
            }
            closed = Some(
                serde_json::json!({"basis":frontier,"target_then_derivative_rows":audits,"conditions":all_conditions}),
            );
            stop = "closed";
            break;
        }
        requested = needed;
        program = combined;
    }
    let result = serde_json::json!({"status":stop,"closed":closed,"scope":"validation-only original-reachable active frontier with persistent exact requested-point history; native rays with point fallback, native proof replay, native exact derivative sums; no alternate elimination or manual master list","input":input_path,"original_targets":targets.iter().map(display).collect::<Vec<_>>(),"physical_input":args[1],"historical_proof_replay":old_proof_replay,"measure_id":sources.measure_id(),"source_count":sources.sources().len(),"shifted_ordinary_slots":shifted,"max_requested":max_requested,"max_rounds":max_rounds,"max_frontier":max_frontier,"max_rules":max_rules,"zero_attempts_per_round":zero_attempts,"depth":depth,"ray_domains_per_requested_point":ray_budget,"point_domains_per_fallback":point_budget,"history":history,"conditions_encountered":all_conditions});
    std::fs::write(
        out.join("result.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    println!("status={stop} rounds={}", history.len());
}
