#!/usr/bin/env python3
"""Exact validation-only two-forest channel evidence from supplied incidence.

This derives a sufficient positive-eta uncut contour certificate. It is not
numerical integration, weighted-source admission, or an eta=0 limit proof.
"""
import hashlib
import json
import copy
from fractions import Fraction
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly/massless-channels.json"


def connected(vertices, edges):
    vertices = set(vertices)
    if not vertices:
        return False
    reached = {min(vertices)}
    while True:
        before = len(reached)
        for tail, head in edges:
            if tail in vertices and head in vertices:
                if tail in reached:
                    reached.add(head)
                if head in reached:
                    reached.add(tail)
        if len(reached) == before:
            return reached == vertices


def certificate(definition, cuts):
    vertices = definition["vertices"]
    assert vertices <= 20, "explicit validation partition budget"
    physical = definition["edges"]
    edges = [edge["vertices"] for i, edge in enumerate(physical) if i not in cuts]
    assert connected(range(vertices), edges)
    orientations = []
    for slot in cuts:
        edge = physical[slot]
        assert Fraction(edge["mass_squared"]) == 0
        potential = sum((charge * Fraction(mu) for charge, mu in zip(edge["charges"], definition["chemical_potentials"])), Fraction(0))
        assert potential != 0
        orientations.append(1 if potential > 0 else -1)
    channels = []
    # Fix vertex zero in the first side to remove complement duplicates.
    for mask in range(1 << (vertices - 1)):
        side = {0} | {v for v in range(1, vertices) if mask & (1 << (v - 1))}
        other = set(range(vertices)) - side
        if not other or not connected(side, edges) or not connected(other, edges):
            continue
        coefficients = [orientation * (int(physical[slot]["vertices"][1] in side) - int(physical[slot]["vertices"][0] in side)) for slot, orientation in zip(cuts, orientations)]
        nonzero = [coefficient for coefficient in coefficients if coefficient]
        if not nonzero:
            kind = "zero"
        elif len(nonzero) == 1:
            kind = "single_lightlike"
        elif len(nonzero) == 2 and sum(nonzero) == 0:
            kind = "difference_of_future_lightlike"
        else:
            kind = "not_certified"
        channels.append({"vertex_partition":[sorted(side),sorted(other)],"future_cut_coefficients":coefficients,"kind":kind})
    return {"cut_slots":cuts,"future_orientations":orientations,"two_forest_partitions":channels,"certified":all(channel["kind"] != "not_certified" for channel in channels)}


def forest_partitions(definition, cuts):
    """Independent finite spanning-forest enumeration checks partition coverage."""
    vertices = definition["vertices"]
    available = [slot for slot in range(len(definition["edges"])) if slot not in cuts]
    partitions = set()
    for chosen in combinations(available, vertices - 2):
        components = [{vertex} for vertex in range(vertices)]
        acyclic = True
        for slot in chosen:
            tail, head = definition["edges"][slot]["vertices"]
            ti = next(i for i, group in enumerate(components) if tail in group)
            hi = next(i for i, group in enumerate(components) if head in group)
            if ti == hi:
                acyclic = False
                break
            merged = components[ti] | components[hi]
            components = [group for i, group in enumerate(components) if i not in (ti, hi)] + [merged]
        if acyclic:
            assert len(components) == 2
            side = next(group for group in components if 0 in group)
            partitions.add(tuple(sorted(side)))
    return partitions


def verify_definition_and_cases(definition, cases):
    vertices, loops = definition["vertices"], definition["loops"]
    assert len(definition["edges"]) - vertices + 1 == loops
    for vertex in range(vertices):
        for loop in range(loops):
            assert sum((Fraction(edge["routing"][loop]) * (int(edge["vertices"][1] == vertex) - int(edge["vertices"][0] == vertex)) for edge in definition["edges"]), Fraction(0)) == 0
    for edge in definition["edges"]:
        for species, charge in enumerate(edge["charges"]):
            assert sum((Fraction(value) * definition["loop_charges"][loop][species] for loop, value in enumerate(edge["routing"])), Fraction(0)) == charge
    for case in cases:
        cuts = case["cut_slots"]
        assert {tuple(channel["vertex_partition"][0]) for channel in case["two_forest_partitions"]} == forest_partitions(definition, cuts)
        # Reverse every edge separately, including its momentum and charge.
        # Future physical channels must remain exactly the same.
        for slot in range(len(definition["edges"])):
            reversed_definition = copy.deepcopy(definition)
            edge = reversed_definition["edges"][slot]
            edge["vertices"].reverse()
            edge["routing"] = [str(-Fraction(value)) for value in edge["routing"]]
            edge["charges"] = [-charge for charge in edge["charges"]]
            reversed_case = certificate(reversed_definition, cuts)
            assert reversed_case["two_forest_partitions"] == case["two_forest_partitions"]


def main():
    manifest = json.loads((ROOT / "docs/finite-density-acceptance.json").read_text())
    records = []
    for family in manifest["four_loop_families"]:
        path = ROOT / family["definition"]
        definition = json.loads(path.read_text())
        # Cut sets come from previously verified incidence/routing certificates;
        # channels themselves use the actual signed incidence, never routing rank.
        cases = [certificate(definition, item["slots"]) for item in family["graph_certificate"]["connected_admissible_occupied_cutsets"]]
        verify_definition_and_cases(definition, cases)
        passed = sum(case["certified"] for case in cases)
        records.append({"family":family["id"],"definition":family["definition"],"definition_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"certified_cutsets":passed,"total_occupied_cutsets":len(cases),"all_cutsets_certified":passed == len(cases),"validation_checks":{"incidence_routing_and_charge_conservation":"passed","explicit_two_forest_coverage":"passed","single_edge_orientation_reversals":len(cases)*len(definition["edges"])},"cases":cases})
        print(f"{family['id']}: {passed}/{len(cases)} occupied cut sets satisfy the sufficient channel certificate")
    assert records[0]["all_cutsets_certified"]
    assert not records[1]["all_cutsets_certified"]
    assert records[2]["all_cutsets_certified"]
    result = {
        "schema":1,
        "status":"exact_graph_channel_evidence_only",
        "certificate":"Every spanning two-forest channel is zero, a single future massless occupied momentum, or a difference of two such momenta; therefore each channel squared is nonpositive in Minkowski signature.",
        "consequence":"For alpha_e>0, sum alpha_e=1, and a common positive uncut auxiliary squared mass eta, the second Symanzik polynomial satisfies F>=eta*U>0. This admits the complete uncut amplitude continuation for eta>0.",
        "not_established":["massless shell/lower-endpoint distribution extension","UV/IR convergence or meromorphic regulator independence","eta=0 physical endpoint or interchange of limits","weighted reduction closure","native numerical acceptance"],
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "families":records,
    }
    REPORT.write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
