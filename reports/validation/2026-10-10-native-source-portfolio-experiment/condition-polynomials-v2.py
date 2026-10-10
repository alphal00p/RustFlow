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
                assert not b or (len(b) == 1 and () in b), "polynomial exponent must be constant"
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

