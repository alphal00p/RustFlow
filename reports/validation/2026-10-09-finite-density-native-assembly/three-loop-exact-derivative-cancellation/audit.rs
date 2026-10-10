use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use symbolica::prelude::*;
use serde_json::{Value, json};

type Label = Vec<i64>;
fn label(v: &Value) -> Label { v.as_array().unwrap().iter().map(|x| x.as_i64().unwrap()).collect() }
fn atom(v: &str) -> Atom { Atom::parse(v, "three_loop_audit", Default::default()).unwrap() }
fn main() {
    let args: Vec<_> = std::env::args().collect();
    let summed_reachability = args.get(4).map(|s| s != "raw").unwrap_or(true);
    let combined_targets = args.get(4).map(|s| s == "targets").unwrap_or(false);
    let saved: Value = serde_json::from_slice(&fs::read(&args[1]).unwrap()).unwrap();
    let roots: Value = serde_json::from_slice(&fs::read(&args[2]).unwrap()).unwrap();
    let rows: BTreeMap<Label, &Value> = saved["reductions"].as_array().unwrap().iter()
        .map(|r| (label(&r["requested"]), r)).collect();
    let global: BTreeSet<Label> = saved["remaining_labels"].as_array().unwrap().iter().map(label).collect();
    let roots: Vec<Label> = roots["requested"].as_array().unwrap().iter().map(label).collect();
    let start = std::time::Instant::now();
    let mut cache: BTreeMap<Label, BTreeMap<Label, Atom>> = BTreeMap::new();
    let mut reduction = |index: &Label| -> Option<BTreeMap<Label, Atom>> {
        if let Some(r) = cache.get(index) { return Some(r.clone()); }
        let out = if let Some(r) = rows.get(index) {
            assert!(r["terminal_terms"].as_array().unwrap().is_empty());
            let mut out = BTreeMap::<Label, Atom>::new();
            for t in r["unresolved"].as_array().unwrap() {
                assert_eq!(t["reason"], "NoApplicableRule");
                *out.entry(label(&t["indices"])).or_insert_with(Atom::new) += atom(t["coefficient"].as_str().unwrap());
            }
            out.retain(|_, c| { *c = c.together().cancel(); !c.is_zero() });
            out
        } else if global.contains(index) {
            [(index.clone(), Atom::num(1))].into_iter().collect()
        } else { return None; };
        cache.insert(index.clone(), out.clone()); Some(out)
    };
    let mut basis = BTreeSet::new();
    for root in &roots { basis.extend(reduction(root).unwrap().into_keys()); }
    let separate_initial = basis.len();
    let mut target_audits = Vec::new();
    if combined_targets {
        assert_eq!(roots.len(), 3);
        basis.clear();
        // All three exact native target coefficients have the same minus Wick
        // phase: five scalar physical powers; or six physical powers and one
        // Gram numerator / two energy numerators. A common minus cannot change
        // the support of an exact linear combination.
        for target_roots in [vec![0usize], vec![1, 2]] {
            let mut sum = BTreeMap::<Label, Atom>::new();
            for i in target_roots.iter().copied() {
                for (leaf, c) in reduction(&roots[i]).unwrap() {
                    *sum.entry(leaf).or_insert_with(Atom::new) -= c;
                }
            }
            let before: BTreeSet<_> = sum.keys().cloned().collect();
            sum.retain(|_, c| { *c = c.together().cancel(); !c.is_zero() });
            let after: BTreeSet<_> = sum.keys().cloned().collect();
            target_audits.push(json!({"original_root_indices":target_roots,"before_sum":before.len(),
                "after_sum":after.len(),"cancelled":before.difference(&after).collect::<Vec<_>>() }));
            basis.extend(after);
        }
    }
    let initial = basis.len();
    let mut done = BTreeSet::new();
    let mut missing = BTreeSet::new();
    let mut changed = Vec::new();
    let mut stages = Vec::new();
    let mut complete_derivatives = 0;
    let mut raw_terms = 0;
    let mut summed_terms = 0;
    for generation in 0..100 {
        let fresh: Vec<_> = basis.difference(&done).cloned().collect();
        if fresh.is_empty() { break; }
        let before = basis.len();
        let mut raw_union = BTreeSet::new();
        let mut summed_union = BTreeSet::new();
        let mut incomplete_rows = 0;
        for parent in fresh {
            done.insert(parent.clone());
            let mut sum = BTreeMap::<Label, Atom>::new();
            let mut components = Vec::new();
            let mut complete = true;
            for slot in [1usize, 2, 4] {
                if parent[slot] == 0 { continue; }
                let mut child = parent.clone(); child[slot] += 1;
                let Some(r) = reduction(&child) else { missing.insert(child); complete = false; continue; };
                for (leaf, coefficient) in r {
                    let coefficient = coefficient * Atom::num(parent[slot]);
                    components.push(json!({"shifted_slot": slot, "leaf": leaf, "coefficient": coefficient.to_string()}));
                    *sum.entry(leaf).or_insert_with(Atom::new) += coefficient;
                }
            }
            if !complete { incomplete_rows += 1; continue; }
            complete_derivatives += 1;
            raw_terms += sum.len();
            let before_sum: BTreeSet<_> = sum.keys().cloned().collect();
            raw_union.extend(before_sum.iter().cloned());
            sum.retain(|_, c| { *c = c.together().cancel(); !c.is_zero() });
            summed_terms += sum.len();
            let after_sum: BTreeSet<_> = sum.keys().cloned().collect();
            let cancelled: Vec<_> = before_sum.difference(&after_sum).cloned().collect();
            if !cancelled.is_empty() {
                changed.push(json!({"generation":generation,"parent":parent,"raw_unique_leaves":before_sum.len(),
                    "summed_unique_leaves":after_sum.len(),"cancelled":cancelled,"components":components}));
            }
            summed_union.extend(after_sum.iter().cloned());
            basis.extend(if summed_reachability { after_sum } else { before_sum });
        }
        stages.push(json!({"generation":generation,"new_basis":basis.len()-before,"basis":basis.len(),
            "raw_derivative_leaf_union":raw_union.len(),"summed_derivative_leaf_union":summed_union.len(),
            "incomplete_derivative_rows":incomplete_rows,"missing_derivative_reductions":missing.len()}));
        if generation == 99 { panic!("generation budget"); }
    }
    let report = json!({"scope":"exact Symbolica together().cancel() within each actual native d_eta integral combination; no new native discovery; incomplete saved derivative rows are skipped, so reachable set is partial and no closure is claimed",
        "summed_reachability":summed_reachability,"combined_original_targets":combined_targets,"target_audits":target_audits,
        "native_program_rules": saved["combined_rules"], "original_labels":roots.len(), "initial_separate_root_leaves":separate_initial,"initial_active_leaves":initial,
        "known_reachable_basis":basis.len(),"historical_leaf_union":global.len(),"complete_derivative_rows":complete_derivatives,
        "incomplete_unique_derivative_labels":missing.len(),"derivative_leaf_occurrences_before_sum":raw_terms,
        "derivative_leaf_occurrences_after_sum":summed_terms,"derivatives_with_exact_cancellation":changed.len(),
        "stages":stages,"changed_derivatives":changed,"missing_labels":missing,
        "conditions":"all original saved native conditions remain required; exact cancellation does not remove their admission obligations",
        "saved_condition_count":saved["conditions"].as_array().map(Vec::len), "seconds":start.elapsed().as_secs_f64()});
    fs::write(&args[3], serde_json::to_vec_pretty(&report).unwrap()).unwrap();
    let mut summary = report.as_object().unwrap().clone(); summary.remove("changed_derivatives");summary.remove("missing_labels");
    println!("{}", serde_json::to_string_pretty(&summary).unwrap());
}
