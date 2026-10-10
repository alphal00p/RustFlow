#!/usr/bin/env python3
"""Prepare the existing bounded native closure pilot with portfolio accounting."""
import hashlib
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
PRIOR = ROOT / "reports/validation/2026-10-10-global-boost-source-experiment/native-active-pilot/pilot.rs"
OUT = BASE / "active-pilot"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert not OUT.exists()
OUT.mkdir()
source = PRIOR.read_text()
assert 'old_proofs_reused' in source
source = source.replace(
    'let mut point_status=None;',
    'let mut point_status=None; let mut portfolio_stats=Vec::new();',
)
needle = 'let exact=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();point_status='
replacement = '''let exact=GuardedProgram::new(sources.clone(),found.rules,[]).unwrap();
    for rule in exact.rules() { if let Some(s)=rule.candidate().stats.guarded_portfolio {
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
    }}
    point_status='''
assert source.count(needle) == 1
source = source.replace(needle, replacement)
needle = '"point_status":point_status,"gaps":gaps'
assert source.count(needle) == 1
source = source.replace(needle, '"point_status":point_status,"portfolio":portfolio_stats,"gaps":gaps')
(OUT / "pilot.rs").write_text(source)
(OUT / "source-binding.json").write_text(json.dumps({
    "scope": "Source-only fresh discovery; same bounded active-frontier algorithm with point portfolio statistics added. No prior rules reused.",
    "base_source": str(PRIOR.relative_to(ROOT)),
    "base_source_sha256": digest(PRIOR),
    "source_sha256": digest(OUT / "pilot.rs"),
    "generator_sha256": digest(Path(__file__)),
}, indent=2) + "\n")
