#!/usr/bin/env python3
"""Definition-only E8 common-energy residue diagnostic; no amplitude values.

Enumerate physical loop bases by their complementary spanning trees, check
exact energy routing and large-one-energy decay of both original numerators.
The occupation-cell observation uses the complete common energy-residue sum,
not independently prescribed cut denominators.
"""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as Q
import json
import hashlib

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "examples/finite_density/five_vertex_eight_edge.json"
OUTPUT = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly/e8-common-energy-basis-diagnostic.json"


def inverse(matrix):
    n = len(matrix)
    rows = [[Q(x) for x in row] + [Q(i == j) for j in range(n)] for i, row in enumerate(matrix)]
    for i in range(n):
        pivot = next((j for j in range(i, n) if rows[j][i]), None)
        if pivot is None:
            return None
        rows[i], rows[pivot] = rows[pivot], rows[i]
        scale = rows[i][i]
        rows[i] = [x / scale for x in rows[i]]
        for j in range(n):
            if i != j:
                factor = rows[j][i]
                rows[j] = [x - factor * y for x, y in zip(rows[j], rows[i])]
    return [row[n:] for row in rows]


def spanning_tree(edges, vertices, slots):
    if len(slots) != vertices - 1:
        return False
    seen = {0}
    while True:
        old = len(seen)
        for slot in slots:
            a, b = edges[slot]["vertices"]
            if a in seen or b in seen:
                seen.update((a, b))
        if len(seen) == old:
            return len(seen) == vertices


def main():
    definition = json.loads(INPUT.read_text())
    edges = definition["edges"]
    assert definition["loops"] == 4 and definition["vertices"] == 5 and len(edges) == 8
    assert definition["targets"] == [
        {"powers": [1] * 8, "numerator": "g1_3^2"},
        {"powers": [2] + [1] * 7, "numerator": "g1_2^2+g1_3*g2_4"},
    ]
    charged = [i for i, e in enumerate(edges) if e["charges"] != [0]]
    neutral = [i for i, e in enumerate(edges) if e["charges"] == [0]]
    assert charged == [0, 2, 4, 6, 7] and neutral == [1, 3, 5]
    assert all(edges[i]["charges"] == [1] for i in charged)
    routing = [[Q(x) for x in e["routing"]] for e in edges]
    bases = []
    for basis in combinations(range(len(edges)), definition["loops"]):
        complement = sorted(set(range(len(edges))) - set(basis))
        inv = inverse([routing[i] for i in basis])
        tree = spanning_tree(edges, definition["vertices"], complement)
        assert (inv is not None) == tree
        if not tree:
            continue
        rows = [[sum(routing[e][k] * inv[k][j] for k in range(4)) for j in range(4)] for e in range(8)]
        basis_charged = sorted(set(basis) & set(charged))
        assert basis_charged
        directions = []
        for j in range(4):
            direction = [inv[k][j] for k in range(4)]
            support = [i for i in range(8) if rows[i][j]]
            # Exact upper polynomial degree: spatial dot products are constants
            # in this energy variable, and P_a0 P_b0 has the displayed degree.
            v = [int(x != 0) for x in direction]
            degree_i91 = 2 * (v[0] + v[2])
            degree_supplemental = max(2 * (v[0] + v[1]), sum(v))
            denominator_degree = 2 * len(support)
            assert denominator_degree - degree_i91 >= 2
            assert denominator_degree - degree_supplemental >= 2
            directions.append({"original_energy_direction": list(map(str, direction)), "physical_edge_support": support, "simple_base_denominator_degree": denominator_degree, "I91_numerator_degree_upper_bound": degree_i91, "supplemental_numerator_degree_upper_bound": degree_supplemental})
        bases.append({"basis_edges": basis, "complementary_spanning_tree": complement, "charged_basis_edges": basis_charged, "directions": directions})
    report = {
        "schema_version": 1,
        "status": "exact_definition_checks_passed",
        "scope": "local occupation-cell cancellation in the complete common energy-residue representation; no four-loop evaluation or contour admission",
        "primary_reference": "https://arxiv.org/html/1609.04339v2#S3",
        "input": str(INPUT.relative_to(ROOT)),
        "input_sha256": hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "physical_basis_count": len(bases),
        "energy_direction_checks": 4 * len(bases),
        "charged_edges": charged,
        "neutral_edges": neutral,
        "occupation_cell": "0<E_e<mu for all five charged edges, with a strict positive margin from each endpoint",
        "cell_argument": "Each L=4 physical energy basis contains a charged edge because only three edges are neutral. Its complete-residue weight contains theta(E_e-mu), which vanishes throughout this open cell. The common energy-integrated sum therefore vanishes there. Raised original-mass differentiation preserves that local zero when performed before specialization, including all grouped terms.",
        "one_energy_arc_checks": "Each original numerator divided by all simple-base denominators decays at least as the inverse square of every physical basis energy with the other basis energies fixed. The actual raised target has stronger denominator decay.",
        "qualification": "This verifies a local consequence of the complete common energy representation and the absence of one-energy polynomial arc obstructions. It does not derive a general thermal distribution for separately evaluated cut amplitudes, establish overlapping endpoint continuation or compare uniform-Feynman and thermal full sums globally.",
        "native_predictions_read": 0,
        "oracle_records_read": 0,
        "bases": bases,
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "physical_basis_count", "energy_direction_checks")}))


if __name__ == "__main__":
    main()
