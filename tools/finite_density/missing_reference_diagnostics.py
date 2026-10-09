#!/usr/bin/env python3
"""Exact definition-only witnesses/kernels for missing references; no amplitudes."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly/missing-reference-diagnostics.json"


def vector(*entries):
    return tuple(Q(entry) for entry in entries)


def add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def scale(factor, value):
    return tuple(factor * component for component in value)


def sub(left, right):
    return add(left, scale(-1, right))


def minkowski(left, right):
    return left[0] * right[0] - sum((a * b for a, b in zip(left[1:], right[1:])), Q(0))


def euclidean(left, right):
    return -minkowski(left, right)


def serialize(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    return value


def e8_witness(definition, q0, q2, q4, q7):
    loops = [q0, sub(q0, q7), q2, sub(q0, q4)]
    edges = [tuple(sum((Q(coefficient) * loops[i][axis] for i, coefficient in enumerate(edge["routing"])), Q(0)) for axis in range(4)) for edge in definition["edges"]]
    charged = [0, 2, 4, 6, 7]
    assert all(edges[slot][0] > 0 and edges[slot][0] < 1 and minkowski(edges[slot], edges[slot]) == 0 for slot in charged)
    neutral = {str(slot): euclidean(edges[slot], edges[slot]) for slot in (1, 3, 5)}
    assert all(value > 0 for value in neutral.values())
    assert add(edges[0], edges[6]) == add(edges[2], edges[4])
    energy_delta = edges[0][0] + edges[6][0] - edges[2][0] - edges[4][0]
    assert energy_delta == 0
    delta_partial_p1_x = edges[0][1] / edges[0][0] - edges[4][1] / edges[4][0]
    assert delta_partial_p1_x == 1
    n91 = euclidean(loops[0], loops[2]) ** 2
    supplemental_base = euclidean(loops[0], loops[1]) ** 2 + euclidean(loops[0], loops[2]) * euclidean(loops[1], loops[3])
    assert n91 != 0 and supplemental_base != 0
    denominator = Q(1)
    for slot in charged:
        denominator *= 2 * edges[slot][0]
    for value in neutral.values():
        denominator *= value
    residues = []
    for missing, sign in [(6, 1), (4, -1), (2, -1), (0, 1)]:
        residues.append({
            "cut_slots": sorted(slot for slot in charged if slot != missing),
            "missing_charged_slot": missing,
            "rho_linear_delta_coefficient": sign * 2 * edges[missing][0],
            "I91_delta_pole_residue": sign * n91 / denominator,
            "supplemental_simple_base_delta_pole_residue": sign * supplemental_base / denominator,
        })
    assert sum(item["I91_delta_pole_residue"] for item in residues) == 0
    assert sum(item["supplemental_simple_base_delta_pole_residue"] for item in residues) == 0
    # This diagnostic uses independent line-pole displacements gamma_e.
    # It is NOT a routed common loop-contour displacement or a thermal proof.
    h_for_equal_line_mass_i0 = sum(Q(sign, 2) / edges[slot][0] for slot, sign in [(0, 1), (6, 1), (2, -1), (4, -1)])
    return {
        "original_minkowski_loop_vectors": loops,
        "routed_edge_vectors": edges,
        "charged_energies": {str(slot): edges[slot][0] for slot in charged},
        "neutral_euclidean_denominators": neutral,
        "I91_numerator": n91,
        "supplemental_simple_base_numerator": supplemental_base,
        "residues": residues,
        "summed_real_delta_pole_residue": "0",
        "delta_partial_original_p1_x": delta_partial_p1_x,
        "uniform_uncut_rho_minus_i0_imaginary_contact_coefficient_in_units_i_pi_delta": {
            "I91": 4 * n91 / denominator,
            "supplemental_simple_base": 4 * supplemental_base / denominator,
        },
        "coherent_line_pole_displacement_H_over_common_mass_i0": h_for_equal_line_mass_i0,
        "raised_target_rule": "Apply -d/dm0_squared to the complete jointly regulated expression with original numerator fixed; include the term where slot0 is uncut. This record does not evaluate the raised amplitude.",
    }


def prism_kernel_check(definition, q1, q2, q3, check_cubic):
    a = minkowski(q1, q1)
    assert minkowski(q2, q2) == 0 and minkowski(q3, q3) == 0
    p1, p2, p4 = sub(q1, q3), q1, sub(q1, q2)
    h12, h13, h23 = euclidean(p4, p4), euclidean(p1, p1), euclidean(sub(p1, p4), sub(p1, p4))
    assert min(h12, h13, h23) > 0
    x0, x1, x2 = Q(1, 2), Q(1, 3), Q(1, 6)
    shifted_k = add(scale(x1, p4), scale(x2, p1))
    loops = [p1, p2, shifted_k, p4]
    edges = [tuple(sum((Q(coefficient) * loops[i][axis] for i, coefficient in enumerate(edge["routing"])), Q(0)) for axis in range(4)) for edge in definition["edges"]]
    assert edges[1] == q1 and edges[5] == q2 and edges[7] == scale(-1, q3)
    assert edges[0] == p1 and edges[3] == p4 and edges[4] == sub(p1, p4)
    assert edges[2] == shifted_k and edges[6] == sub(shifted_k, p4) and edges[8] == sub(p1, shifted_k)
    direct = euclidean(p1, p2) ** 2 + euclidean(p1, shifted_k) * euclidean(p2, p4)
    reduced = (h13 - a) ** 2 / 4 + (h12 - a) / 2 * (x1 * (h13 + h12 - h23) / 2 + x2 * h13)
    assert direct == reduced
    cubic = None
    if check_cubic:
        assert a == 0
        # Gaussian covariance terms vanish because external p2²=0.
        cubic = euclidean(p2, shifted_k) ** 3
        assert cubic == (x1 * h12 + x2 * h13) ** 3 / 8
    return {"cut_mass_squared_a": a,"future_cut_vectors": [q1,q2,q3],"h12_h13_h23": [h12,h13,h23],"feynman_parameters": [x0,x1,x2],"triangle_F": x0*x1*h12+x0*x2*h13+x1*x2*h23,"supplemental_simple_base_direct": direct,"supplemental_simple_base_parameter_kernel": reduced,"I115_parameter_kernel_if_massless": cubic}


def prism_two_cut_parameters():
    checks = []
    for parameters in [(1,2,3,4,5,6),(Q(1,2),Q(1,3),Q(2,5),Q(3,7),Q(5,11),Q(7,13))]:
        a,b,c,d,e,f = map(Q,parameters)
        A,B = a+c+e+f,b+d+f
        U = A*B-f*f
        phi = a*b*c+a*c*d+a*c*f+b*c*f+a*b*d+b*c*d+b*d*e+b*d*f+a*d*f
        psi = e*(a*b+a*d+a*f+b*f)
        for virtuality in (Q(0),Q(-1,4)):
            h = Q(2)
            q_dot_p = (h+virtuality)/2
            b1_squared = c*c*h+e*e*virtuality+2*c*e*q_dot_p
            b2_squared = d*d*h
            b1_dot_b2 = c*d*h+e*d*q_dot_p
            direct = U*((c+d)*h+e*virtuality) - B*b1_squared-2*f*b1_dot_b2-A*b2_squared
            assert direct == h*phi+virtuality*psi
            assert U > 0 and phi > 0 and psi > 0
            checks.append({"parameters":parameters,"external_h":h,"external_q_squared":virtuality,"U":U,"Phi":phi,"Psi":psi,"F_from_Gaussian_completion":direct,"F_from_positive_polynomials":h*phi+virtuality*psi})
    return checks


def main():
    paths = [ROOT / "examples/finite_density/five_vertex_eight_edge.json", ROOT / "examples/finite_density/triangular_prism.json"]
    e8, prism = [json.loads(path.read_text()) for path in paths]
    assert e8["targets"][1] == {"powers": [2,1,1,1,1,1,1,1],"numerator":"g1_2^2+g1_3*g2_4"}
    assert prism["targets"][1] == {"powers": [1,2,1,1,1,1,1,1,1],"numerator":"g1_2^2+g1_3*g2_4"}
    output = {
        "schema":1,
        "status":"exact_definition_diagnostics_passed_no_references_generated",
        "oracle_values_read":0,"amf_predictions_read":0,"generated_amplitude_reference_values":0,
        "input_sha256": {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "e8_interior_pole_witnesses":[
            e8_witness(e8,vector("1/2","1/2",0,0),vector("1/2",0,0,"1/2"),vector("1/2",0,0,"-1/2"),vector("1/4",0,"1/4",0)),
            e8_witness(e8,vector("1/3","1/3",0,0),vector("1/2",0,0,"1/2"),vector("1/4",0,0,"-1/4"),vector("1/5",0,"1/5",0)),
        ],
        "prism_triangle_numerator_checks":[
            prism_kernel_check(prism,vector("1/2","1/2",0,0),vector("1/3",0,"1/3",0),vector("1/4",0,0,"1/4"),True),
            prism_kernel_check(prism,vector("3/4","1/2",0,0),vector("1/3",0,"1/3",0),vector("1/4",0,0,"1/4"),False),
        ],
        "prism_two_cut_parameter_checks":prism_two_cut_parameters(),
        "scope":"Exact interior-pole/residue and original-numerator parameter identities only. Common thermal prescription, analytic subtraction, integrated Laurent values and numerical reference refinements remain required.",
    }
    REPORT.write_text(json.dumps(serialize(output),indent=2)+"\n")
    print("Passed two exact E8 pole/contact witnesses and massless/massive prism numerator kernel checks. No amplitude reference generated.")


if __name__ == "__main__":
    main()
