"""Independent exact Gaussian determinant check of native full-kinematic evidence.

Usage: python compare_full_kinematics.py NATIVE_DIRECTORY [OUTPUT_JSON]
Rational polynomial arithmetic only; no HEPKit, production evaluator, reference
period, SymPy, or supplied numerical answer is imported. Each native routed row
is treated as exact input; compare_native.py separately binds the two-cut rows
to the independent original-routing calculation.
"""
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import gzip
import hashlib
import json
import sys


def add(*polys):
    out = {}
    for poly in polys:
        for e, c in poly.items():
            out[e] = out.get(e, Q(0)) + c
    return {e: c for e, c in out.items() if c}


def scale(poly, coefficient):
    return {e: c * coefficient for e, c in poly.items() if c * coefficient}


def mul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            exponent = tuple(x + y for x, y in zip(e, f))
            out[exponent] = out.get(exponent, Q(0)) + c * d
    return {e: c for e, c in out.items() if c}


def determinant(matrix, n):
    if not matrix:
        return {(0,) * n: Q(1)}
    out = {}
    for order in permutations(range(len(matrix))):
        sign = (-1) ** sum(order[i] > order[j] for i in range(len(order)) for j in range(i + 1, len(order)))
        term = {(0,) * n: Q(sign)}
        for i, j in enumerate(order):
            term = mul(term, matrix[i][j])
        out = add(out, term)
    return out


def decoded(terms):
    return {tuple(t["exponents"]): Q(t["coefficient"]) for t in terms}


def check(support, k, h):
    active = support["active_slots"]
    n = len(active)
    rows = [[Q(a) for a in support["routed_rows"][slot]] for slot in active]
    alpha = [{tuple(int(i == j) for j in range(n)): Q(1)} for i in range(n)]
    moment = lambda i, j: add(*(scale(a, row[i] * row[j]) for a, row in zip(alpha, rows)))
    a = [[moment(k + i, k + j) for j in range(h)] for i in range(h)]
    u = determinant(a, n)
    full = support["full_kinematics"]
    if not u:
        assert full is None, "rank-deficient support carries full-rank polynomials"
        return "rank_deficient"
    assert full is not None
    assert (full["compact_loops"], full["virtual_loops"]) == (k, h)
    assert full["parameter_slots"] == active
    assert decoded(full["u"]) == u, ("U", active)
    adj = [[scale(determinant([[a[r][c] for c in range(h) if c != i]
                              for r in range(h) if r != j], n), (-1) ** (i + j))
            for j in range(h)] for i in range(h)]
    b = [[moment(i, k + j) for j in range(h)] for i in range(k)]
    def external(i, j):
        contracted = add(*(mul(mul(b[i][r], adj[r][s]), b[j][s])
                           for r in range(h) for s in range(h)))
        return add(contracted, scale(mul(u, moment(i, j)), Q(-1)))
    assert len(full["pairs"]) == k * (k - 1) // 2
    for pair in full["pairs"]:
        i, j = pair["pair"]
        expected = external(i, j)
        assert all(c > 0 for c in expected.values())
        assert decoded(pair["terms"]) == expected, ("pair", active, i, j)
    assert len(full["diagonal_masses"]) == k
    for i, terms in enumerate(full["diagonal_masses"]):
        assert decoded(terms) == external(i, i), ("off-null diagonal", active, i)
    assert len(full["virtual_masses"]) == n
    for i, terms in enumerate(full["virtual_masses"]):
        assert decoded(terms) == mul(u, alpha[i]), ("independent virtual mass", active, i)
    return "full_rank"


directory = Path(sys.argv[1])
runs = []
for name, k, h in [(f"prism-cuts-{a}-{b}-supports.json", 2, 2)
                   for a, b in [(1, 5), (1, 7), (5, 7)]] + [
                       ("prism-cuts-1-5-7-one-virtual-supports.json", 3, 1)]:
    path = directory / name
    if not path.exists():
        path = path.with_suffix(".json.gz")
    raw = path.read_bytes()
    payload = gzip.decompress(raw) if path.suffix == ".gz" else raw
    data = json.loads(payload)
    counts = {"rank_deficient": 0, "full_rank": 0}
    for support in data["supports"]:
        counts[check(support, k, h)] += 1
    runs.append({"file": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
                 "supports": len(data["supports"]), **counts})
result = {"scope": "independent exact routed Gaussian determinants versus native full off-null HEPKit U/F; query only; no admission or integral values",
          "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          "passed": sum(r["supports"] for r in runs), "failed": 0, "runs": runs}
output = Path(sys.argv[2]) if len(sys.argv) > 2 else directory / "independent-full-kinematic-comparison.json"
output.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({"passed": result["passed"], "failed": 0, "output": str(output)}))
