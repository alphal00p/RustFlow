#!/usr/bin/env python3
"""Exact rational-associate check of archived polynomial condition strings.

Independent finite sparse-polynomial arithmetic; no eval, floating point,
factorization, root-set inference, or native rule mutation.
"""
import ast
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path


def add(a, b):
    out = dict(a)
    for monomial, value in b.items():
        out[monomial] = out.get(monomial, Q(0)) + value
    return {m: c for m, c in out.items() if c}


def mul(a, b):
    out = {}
    for am, ac in a.items():
        for bm, bc in b.items():
            powers = dict(am)
            for name, power in bm:
                powers[name] = powers.get(name, 0) + power
            monomial = tuple(sorted(powers.items()))
            out[monomial] = out.get(monomial, Q(0)) + ac * bc
    assert len(out) <= 1000
    return {m: c for m, c in out.items() if c}


def parse(expression):
    def visit(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return {} if node.value == 0 else {(): Q(node.value)}
        if isinstance(node, ast.Name):
            return {((node.id, 1),): Q(1)}
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return {m: -c for m, c in visit(node.operand).items()}
            if isinstance(node.op, ast.UAdd):
                return visit(node.operand)
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add):
                return add(a, b)
            if isinstance(node.op, ast.Sub):
                return add(a, {m: -c for m, c in b.items()})
            if isinstance(node.op, ast.Mult):
                return mul(a, b)
            if isinstance(node.op, ast.Div):
                assert len(b) == 1 and () in b and b[()]
                return {m: c / b[()] for m, c in a.items()}
            if isinstance(node.op, ast.Pow):
                exponent = Q(0) if not b else b.get(())
                assert exponent is not None and exponent.denominator == 1 and 0 <= exponent <= 16
                out = {(): Q(1)}
                for _ in range(int(exponent)):
                    out = mul(out, a)
                return out
        raise ValueError(ast.dump(node))
    return visit(ast.parse(expression.replace('^', '**'), mode='eval').body)


def canonical(polynomial):
    if not polynomial:
        return None
    leading = polynomial[min(polynomial)]
    return tuple((m, c / leading) for m, c in sorted(polynomial.items()))


def condition_gate(baseline, trial):
    base = [canonical(parse(s)) for s in baseline]
    matches = []
    for text in trial:
        polynomial = parse(text)
        if polynomial and set(polynomial) == {()}:
            matches.append({"condition": text, "reason": "nonzero rational constant"})
            continue
        normalized = canonical(polynomial)
        match = next((i for i, value in enumerate(base) if value is not None and value == normalized), None)
        if match is None:
            return False, matches + [{"condition": text, "reason": "not rational-associate of baseline condition"}]
        matches.append({"condition": text, "reason": "rational associate", "baseline_condition": baseline[match]})
    return True, matches


assert condition_gate(["eta"], ["-2*eta", "3/5"])[0]
assert condition_gate(["eta"], [])[0]
assert not condition_gate(["eta"], ["eta^2"])[0]  # Deliberately no root-set theorem.
assert not condition_gate(["eta", "epsilon"], ["eta*epsilon"])[0]
assert not condition_gate(["eta"], ["0"])[0]
assert not condition_gate(["eta"], ["epsilon"])[0]

here = Path(__file__).resolve().parent
source = here.parent / 'selector-controls-comparison.json'
data = json.loads(source.read_text())
baseline = data['selectors']['original']['rows']
arms = [data['selectors'][name]['rows'] for name in ('strict', 'index-dependent')]
results = []
for number, row in enumerate(baseline):
    chosen_terms = row['baseline_rhs_terms']
    chosen_arm = 'baseline'
    decisions = []
    for name, arm in zip(('A', 'C'), arms):
        trial = arm[number]
        assert trial['target'] == row['target'] and trial['baseline_conditions'] == row['baseline_conditions']
        passed, matches = condition_gate(row['baseline_conditions'], trial['trial_conditions'])
        shorter = trial['trial_applied'] and trial['trial_rhs_terms'] < chosen_terms
        selected = shorter and passed
        if trial['strictly_shorter']:
            decisions.append({'arm': name, 'rhs_terms': trial['trial_rhs_terms'], 'condition_gate': passed, 'matches': matches, 'selected_at_this_step': selected})
        if selected:
            chosen_terms, chosen_arm = trial['trial_rhs_terms'], name
    results.append({'target': row['target'], 'baseline_rhs_terms': row['baseline_rhs_terms'], 'original_baseline_conditions': row['baseline_conditions'], 'shorter_trial_decisions': decisions, 'selected_arm': chosen_arm, 'selected_rhs_terms': chosen_terms})

total = sum(r['selected_rhs_terms'] for r in results)
assert total == 181
assert sum(r['selected_arm'] != 'baseline' for r in results) == 2
out = {
    'scope': 'Independent exact rational-polynomial condition audit of saved native point rules; no new native search or closure claim',
    'inputs_sha256': {source.name: hashlib.sha256(source.read_bytes()).hexdigest(), Path(__file__).name: hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    'baseline_rhs_terms': 192,
    'two_arm_condition_safe_rhs_terms': total,
    'improved_points': 2,
    'simple_arithmetic_guard_checks': 6,
    'points': results,
}
(here / 'condition-check.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({k: out[k] for k in ('baseline_rhs_terms', 'two_arm_condition_safe_rhs_terms', 'improved_points')}))
