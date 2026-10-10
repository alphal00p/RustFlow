#![allow(dead_code)]
//! Diagnostic only: reconstruct the complete historical source corpus embedded
//! in a locally generated program, then ask native RustRed to replay and search.
use bincode::Decode;
use rustred::persistence::{
    BinaryIoLimits, CoefficientId, DecodedCoefficientTable, SectionTag, inspect_program,
};
use rustred::solver::guarded::{
    GuardedProgram, GuardedSource, GuardedSourceSystem, IndexBounds, IndexDomain, IndexRole,
};
use rustred::solver::{Integral, Power, SearchOptions, Term};
use std::sync::Arc;
type Domain = Vec<(Option<i64>, Option<i64>)>;
type Label = Vec<(bool, i16)>;
#[derive(Decode)]
struct Tr {
    integral: Label,
    coefficient: usize,
}
#[derive(Decode)]
struct Sr {
    id: String,
    domain: Domain,
    conditions: Vec<usize>,
    terms: Vec<Tr>,
}
#[derive(Decode)]
struct Sd {
    row: usize,
    integral: Label,
    shifts: Vec<i16>,
}
#[derive(Decode)]
struct Rr {
    fixed: Vec<Option<i16>>,
    target: Label,
    rhs: Vec<Tr>,
    sources: Vec<Sd>,
    domain: Domain,
    discovery_domain: Domain,
    conditions: Vec<usize>,
    sector: Vec<bool>,
    permutation: Option<Vec<usize>>,
}
#[derive(Decode)]
struct Rec {
    schema: String,
    measure: String,
    roles: Vec<u8>,
    indices: Vec<usize>,
    sources: Vec<Sr>,
    zero_domains: Vec<Domain>,
    rules: Vec<Rr>,
    terminals: Vec<Vec<i64>>,
}
fn domain<const N: usize>(v: &Domain) -> IndexDomain<N> {
    IndexDomain::new(
        v.iter()
            .map(|&(lo, hi)| IndexBounds::new(lo, hi).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
    .unwrap()
}
fn label<const N: usize>(v: &Label) -> Integral<N> {
    Integral::new(
        v.iter()
            .map(|&(s, p)| Power::new(s, p).unwrap())
            .collect::<Vec<_>>()
            .try_into()
            .unwrap(),
    )
}
fn bounds<const N: usize>(v: &IndexDomain<N>) -> Domain {
    v.bounds().iter().map(|b| (b.lower(), b.upper())).collect()
}

use std::collections::{BTreeMap, BTreeSet};
use symbolica::prelude::*;

type Point = [i64; 12];
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
fn reduce(program: &GuardedProgram<12>, row: &Row, admitted: &IndexDomain<12>) -> Reduced {
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
    std::fs::create_dir_all(&out).unwrap();
    let bytes = std::fs::read(&args[1]).unwrap();
    let limits = BinaryIoLimits::default();
    let envelope = inspect_program(&bytes, limits).unwrap();
    let (record, used): (Rec, usize) = bincode::decode_from_slice(
        envelope.section(SectionTag::PROGRAM).unwrap(),
        bincode::config::standard(),
    )
    .unwrap();
    assert_eq!(used, envelope.section(SectionTag::PROGRAM).unwrap().len());
    assert_eq!(record.schema, "rustred.guarded-source-program.v2");
    let table = DecodedCoefficientTable::import_generated_normalized(
        envelope.section(SectionTag::SYMBOLICA_STATE).unwrap(),
        envelope.section(SectionTag::COEFFICIENTS).unwrap(),
        limits,
    )
    .unwrap();
    let polynomial = |id| {
        let c = table
            .coefficient(CoefficientId::try_from_index(id).unwrap())
            .unwrap()
            .clone();
        assert!(c.denominator.is_one());
        c.numerator
    };
    let roles: [IndexRole; 12] = record
        .roles
        .iter()
        .map(|r| match r {
            0 => IndexRole::Ordinary,
            1 => IndexRole::RequiredCut,
            2 => IndexRole::Occupation,
            _ => panic!(),
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    let rows = record
        .sources
        .iter()
        .map(|s| {
            GuardedSource::new(
                s.id.clone(),
                s.terms
                    .iter()
                    .map(|t| Term {
                        integral: label(&t.integral),
                        coefficient: polynomial(t.coefficient),
                    })
                    .collect(),
                domain(&s.domain),
            )
            .with_nonzero_conditions(s.conditions.iter().map(|&id| polynomial(id)).collect())
        })
        .collect();
    let sources = Arc::new(
        GuardedSourceSystem::new(
            record.measure.clone(),
            roles,
            record.indices.clone().try_into().unwrap(),
            rows,
        )
        .unwrap()
        .with_zero_domains(record.zero_domains.iter().map(domain).collect())
        .unwrap(),
    );
    // Old proof containers are deliberately not reused. Fresh discovery binds
    // every rule to exactly the reconstructed physical source corpus. Native v2
    // rejects the old v1 proof semantics after the zero-projection owner changed.
    assert_ne!(std::env::var("PILOT_REPLAY_OLD").as_deref(), Ok("true"));
    let old_proof_replay: Option<usize> = None;

    // Extract the actual separate target sums from the sealed input/family
    // metadata, rather than guessing relative phases or loading any period.
    let input_path = "sealed source identity: original massless_three_loop_chain input";
    let mut extracted: BTreeMap<usize, Row> = BTreeMap::new();
    for text in record.measure.split("target=").skip(1) {
        let Some((ordinal, tail)) = text.split_once(':') else {
            continue;
        };
        let Ok(ordinal) = ordinal.parse::<usize>() else {
            continue;
        };
        let Some(end) = tail.find(']') else { continue };
        let raw: Vec<i64> = serde_json::from_str(&tail[..=end]).unwrap();
        let coefficient = tail[end + 1..]
            .strip_prefix(':')
            .unwrap()
            .chars()
            .take_while(|c| c.is_ascii_digit() || matches!(c, '-' | '+' | '/' | '.'))
            .collect::<String>();
        assert!(!coefficient.is_empty());
        let c = Atom::parse(&coefficient, "point_pilot", Default::default()).unwrap();
        let mut point = [0; 12];
        point[..raw.len()].copy_from_slice(&raw);
        let entry = extracted.entry(ordinal).or_default();
        if let Some(old) = entry.get(&point) {
            assert_eq!(old, &c);
        } else {
            entry.insert(point, c);
        }
    }
    assert_eq!(extracted.len(), 2);
    assert_eq!(extracted.values().map(|r| r.len()).sum::<usize>(), 3);
    let targets = extracted.into_values().collect::<Vec<_>>();
    let geometry_marker = record.measure.split("physical=").nth(1).unwrap();
    let physical_slots: usize = geometry_marker.split(';').next().unwrap().parse().unwrap();
    let input_slots = roles
        .iter()
        .position(|&r| r == IndexRole::Occupation)
        .unwrap();
    let physical_arity = roles
        .iter()
        .rposition(|&r| r != IndexRole::Ordinary)
        .unwrap()
        + 1;
    let shifted = (0..physical_slots)
        .filter(|&s| roles[s] == IndexRole::Ordinary)
        .collect::<Vec<_>>();
    assert_eq!(shifted, vec![1, 2, 3, 4]);
    let mut b = *IndexDomain::for_roles(&roles).bounds();
    for x in &mut b[physical_slots..input_slots] {
        *x = IndexBounds::new(None, Some(0)).unwrap();
    }
    b[physical_arity..].fill(IndexBounds::fixed(0));
    let admitted = IndexDomain::new(b).unwrap();
    let max_requested = setting("PILOT_REQUESTS", 4096);
    let max_rounds = setting("PILOT_ROUNDS", 8);
    let max_frontier = setting("PILOT_FRONTIER", 128);
    let max_rules = setting("PILOT_RULES", 32768);
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
    let configuration = serde_json::json!({"request_history_policy":"explicit requested labels only; basis and raw derivative labels deferred until requested; max_requested caps cumulative history","source_program":args[1],"source_count":sources.sources().len(),"original_targets":targets.iter().map(display).collect::<Vec<_>>(),"old_proofs_reused":false,"old_proof_optional_replay":old_proof_replay,"max_requested":max_requested,"max_rounds":max_rounds,"max_frontier":max_frontier,"max_rules":max_rules,"zero_attempts_per_round":zero_attempts,"depth":depth,"ray_domains_per_requested_point":ray_budget,"point_domains_per_fallback":point_budget,"shifted_ordinary_slots":shifted,"measure_id":sources.measure_id()});
    std::fs::write(
        out.join("configuration.json"),
        serde_json::to_vec_pretty(&configuration).unwrap(),
    )
    .unwrap();
    let mut history = Vec::<serde_json::Value>::new();
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
    let mut first_round = 0;
    let mut resume_evidence = serde_json::Value::Null;
    if let Some(parent_arg) = args.get(3) {
        let parent = std::path::Path::new(parent_arg);
        let read = |name: &str| -> serde_json::Value {
            serde_json::from_slice(&std::fs::read(parent.join(name)).unwrap()).unwrap()
        };
        let resource = read("resources.json");
        let run_binding = read("run-binding.json");
        assert_eq!(resource["exit_code"], 124);
        assert_eq!(run_binding["exit_code"], 124);
        assert_eq!(run_binding["post_run_inputs_unchanged"], true);
        assert_eq!(
            run_binding["engine"],
            if cfg!(source_portfolio) {
                "proposal"
            } else {
                "baseline"
            }
        );
        assert!(
            !parent.join("result.json").exists() && !parent.join("closed.bin").exists(),
            "only unfinished timeout checkpoints may resume"
        );
        assert_eq!(
            read("configuration.json"),
            configuration,
            "source identity, targets and every limit must match"
        );
        assert_eq!(
            (max_rounds, max_requested, max_frontier, max_rules),
            (16, 4096, 1024, 65536)
        );
        assert_eq!(
            (zero_attempts, depth, ray_budget, point_budget),
            (8192, 3, 1, 1)
        );
        history = serde_json::from_value(read("progress.json")).unwrap();
        assert!(!history.is_empty());
        for (i, item) in history.iter().enumerate() {
            assert_eq!(item["round"].as_u64(), Some(i as u64));
        }
        let last = history.last().unwrap();
        let last_round = last["round"].as_u64().unwrap() as usize;
        first_round = last_round + 1;
        assert!(
            first_round < max_rounds,
            "global round budget already exhausted"
        );
        let round_file = format!("round-{last_round:03}.json");
        let proof_file = format!("round-{last_round:03}.bin");
        let state = read(&round_file);
        assert_eq!(state["round"].as_u64(), Some(last_round as u64));
        assert_eq!(state["frontier_budget_exhausted"], false);
        let points = |field: &str| -> BTreeSet<Point> {
            let raw: Vec<Point> = serde_json::from_value(state[field].clone()).unwrap();
            let set = raw.iter().copied().collect::<BTreeSet<_>>();
            assert_eq!(set.len(), raw.len());
            assert!(set.iter().all(|i| admitted.contains(i)));
            set
        };
        frontier = points("frontier");
        attempted_points = points("attempted_exact_points");
        requested = points("next_needed");
        assert!(roots.is_subset(&requested) && frontier.is_subset(&requested));
        assert!(points("requested").is_subset(&attempted_points));
        assert!(
            !points("pending_new_derivative_labels").is_empty(),
            "a completed closing round must not resume"
        );
        assert!(
            frontier.len() <= max_frontier
                && attempted_points.union(&requested).count() <= max_requested
        );
        all_conditions = serde_json::from_value(state["conditions_encountered"].clone()).unwrap();
        program = GuardedProgram::decode_generated(
            &std::fs::read(parent.join(&proof_file)).unwrap(),
            sources.clone(),
            limits,
        )
        .unwrap();
        assert!(
            program.terminals().is_empty(),
            "a resumable program cannot contain declared master terminals"
        );
        assert!(program.rules().len() <= max_rules);
        assert_eq!(
            state["rule_count"].as_u64(),
            Some(program.rules().len() as u64)
        );
        assert_eq!(last["rules"], state["rule_count"]);
        assert_eq!(last["frontier"].as_u64(), Some(frontier.len() as u64));
        assert_eq!(
            last["requested"].as_u64(),
            Some(points("requested").len() as u64)
        );
        let mut replayed_targets = Vec::new();
        for row in &targets {
            let replay = reduce(&program, row, &admitted);
            assert!(replay.failures.is_empty());
            assert!(replay.conditions.is_subset(&all_conditions));
            replayed_targets.push(display(&replay.leaves));
        }
        assert_eq!(
            serde_json::json!(replayed_targets),
            state["original_target_rows"]
        );
        for row in state["derivative_rows"].as_array().unwrap() {
            let basis: Point = serde_json::from_value(row["basis"].clone()).unwrap();
            assert!(frontier.contains(&basis));
            let raw = derivative(basis, &shifted);
            assert!(raw.keys().all(|i| attempted_points.contains(i)));
            assert_eq!(serde_json::json!(display(&raw)), row["raw_derivative"]);
            let replay = reduce(&program, &raw, &admitted);
            assert!(replay.failures.is_empty());
            assert!(replay.conditions.is_subset(&all_conditions));
            assert_eq!(
                serde_json::json!(display(&replay.leaves)),
                row["reduced_actual_derivative"]
            );
        }
        resume_evidence = serde_json::json!({"parent":parent_arg,"last_completed_round":last_round,
   "first_continuation_round":first_round,"restored_rules":program.rules().len(),
   "restored_frontier":frontier.len(),"restored_attempted_exact_points":attempted_points.len(),
   "restored_requested":requested.len(),"restored_conditions":all_conditions,
   "native_decode_and_saved_target_derivative_replay":true,
   "original_corpus_rules_imported":false,"only_this_engine_parent_discovered_rules_restored":true});
        std::fs::write(
            out.join("resume-state.json"),
            serde_json::to_vec_pretty(&resume_evidence).unwrap(),
        )
        .unwrap();
        std::fs::write(
            out.join("progress.json"),
            serde_json::to_vec_pretty(&history).unwrap(),
        )
        .unwrap();
    }
    for round in first_round..max_rounds {
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
        let mut combined = zero_program.union_replayed(program, [], max_rules).unwrap();
        let mut fresh_programs = Vec::new();
        let mut probes = Vec::new();
        for point in points {
            let existing = combined.apply(&point).unwrap();
            if !matches!(
                existing.status,
                rustred::solver::guarded::GuardedApplicationStatus::Unresolved(_)
            ) {
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
            let fallback = matches!(
                ray_status,
                rustred::solver::guarded::GuardedApplicationStatus::Unresolved(_)
            );
            fresh_programs.push(ray_program);
            let mut point_status = None;
            let mut portfolio_stats = Vec::<serde_json::Value>::new();
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
                #[cfg(source_portfolio)]
                for rule in exact.rules() {
                    if let Some(s) = rule.candidate().stats.guarded_portfolio {
                        portfolio_stats.push(serde_json::json!({
      "baseline_rhs_terms":s.baseline_rhs_terms,"selected_rhs_terms":s.selected_rhs_terms,
      "selected_arm":s.selected_arm,"completion":s.completion,
      "baseline_elapsed":s.baseline_elapsed.as_secs_f64(),
      "policy_elapsed":s.policy_elapsed.as_secs_f64(),
      "preseal_elapsed":s.preseal_elapsed.as_secs_f64(),
      "trials":s.trials.iter().map(|t|serde_json::json!({
       "selector":t.selector,"attempted_rows":t.attempted_rows,"accepted_rows":t.accepted_rows,
       "guard_rejected_rows":t.guard_rejected_rows,"empty_rows":t.empty_rows,
       "rhs_terms":t.rhs_terms,"completion":t.completion,
       "exact_trace_rows":t.exact_trace_rows,"exact_trace_terms":t.exact_trace_terms,
       "search_elapsed":t.search_elapsed.as_secs_f64()
      })).collect::<Vec<_>>()
     }));
                    }
                }
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
                    left.union_replayed(right, [], max_rules).unwrap()
                } else {
                    left
                });
            }
            fresh_programs = next;
        }
        if let Some(fresh) = fresh_programs.pop() {
            combined = combined.union_replayed(fresh, [], max_rules).unwrap();
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
            std::fs::write(
                out.join("closed.bin"),
                audited.encode_native(limits).unwrap(),
            )
            .unwrap();
            closed = Some(
                serde_json::json!({"basis":frontier,"target_then_derivative_rows":audits,"conditions":all_conditions}),
            );
            stop = "closed";
            break;
        }
        requested = needed;
        program = combined;
    }
    let result = serde_json::json!({"status":stop,"closed":closed,"resume":resume_evidence,"scope":"validation-only original-reachable active frontier with persistent exact requested-point history; native rays with point fallback, native proof replay, native exact derivative sums; no alternate elimination or manual master list","input":input_path,"original_targets":targets.iter().map(display).collect::<Vec<_>>(),"historical_source_program":args[1],"historical_proof_replay":old_proof_replay,"measure_id":sources.measure_id(),"source_count":sources.sources().len(),"shifted_ordinary_slots":shifted,"max_requested":max_requested,"max_rounds":max_rounds,"max_frontier":max_frontier,"max_rules":max_rules,"zero_attempts_per_round":zero_attempts,"depth":depth,"ray_domains_per_requested_point":ray_budget,"point_domains_per_fallback":point_budget,"history":history,"conditions_encountered":all_conditions});
    std::fs::write(
        out.join("result.json"),
        serde_json::to_vec_pretty(&result).unwrap(),
    )
    .unwrap();
    println!("status={stop} rounds={}", history.len());
}
