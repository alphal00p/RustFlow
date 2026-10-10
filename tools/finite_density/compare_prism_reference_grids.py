#!/usr/bin/env python3
"""Validate independent prism reference refinements from saved samples only.

No AMF algorithm or prediction is imported. The published uncertainty is an
empirical envelope, not a rigorous quadrature or interpolation error bound.
"""
import argparse
import gzip
import hashlib
import json
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROFILES = (
    ("laurent-baseline", 70, 10000, "1/30", 14),
    ("laurent-nearby-grid", 70, 500000, "1/30", 14),
    ("laurent-fine-baseline", 70, 1000000, "1/30", 14),
    ("laurent-fine-precision", 80, 1000000, "1/30", 14),
    ("laurent-fine-quadrature", 80, 1000000, "1/40", 18),
)
FIELDS = ("full_raw_value", "full_msbar_reference_normalization")
ORDERS = range(-4, 1)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def relative_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def complex_value(text):
    require(text.startswith("(") and text.endswith("i)"), "complex value format changed")
    real, imag = text[1:-2].split("+", 1)
    return Decimal(real), Decimal(imag)


def interpolate(samples, field, count=12):
    """Fit epsilon^6 I with exact rational Lagrange weights, without forced zeros."""
    points = [Fraction(row["epsilon"]) for row in samples[:count]]
    coefficients = [[Decimal(0), Decimal(0)] for _ in points]
    for x, row in zip(points, samples):
        poly, denominator = [Fraction(1)], Fraction(1)
        for other in points:
            if other == x:
                continue
            product = [Fraction(0)] * (len(poly) + 1)
            for k, coefficient in enumerate(poly):
                product[k] -= other * coefficient
                product[k + 1] += coefficient
            poly = product
            denominator *= x - other
        value = complex_value(row[field])
        for k, weight in enumerate(poly):
            weight = weight / denominator * x**6
            weight = Decimal(weight.numerator) / Decimal(weight.denominator)
            for part in range(2):
                coefficients[k][part] += weight * value[part]
    return {str(k - 6): {"re": str(v[0]), "im": str(v[1])} for k, v in enumerate(coefficients)}


def compare(a, b):
    distance = max(abs(Decimal(a[k]) - Decimal(b[k])) for k in ("re", "im"))
    scale = max(abs(Decimal(v[k])) for v in (a, b) for k in ("re", "im"))
    small = scale < Decimal("1e-20")
    tolerance = Decimal("1e-25") if small else scale * Decimal("1e-12")
    return {"absolute_difference": str(distance), "relative_difference": str(distance / scale) if scale else "0", "criterion": "absolute" if small else "relative", "tolerance": str(tolerance), "passed": distance <= tolerance}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-incomplete", action="store_true")
    parser.add_argument("--reference-output", type=Path, help="write the independent reference only after every gate passes")
    args = parser.parse_args()
    directory, output = args.directory.resolve(), args.output.resolve()
    files = set()
    reference_sources = None
    executable = None
    profiles, pending = [], []
    checks = {"analytic_regulator_residue_cancellations": 0, "sample_contour_pole_checks": 0, "uniform_stirling_tail_checks": 0, "archived_contours_replayed": 0}
    fixture = ROOT / "examples/finite_density/triangular_prism.json"
    definition = read(fixture)
    require(definition["loops"] == 4 and definition["vertices"] == 6, "prism graph dimension changed")
    require(definition["loop_charges"] == [[0], [1], [0], [0]], "prism loop chemical assignments changed")
    require([e["charges"] for e in definition["edges"]] == [[0], [1], [0], [0], [0], [1], [0], [-1], [0]], "prism line chemical assignments changed")
    require([e["vertices"] for e in definition["edges"]] == [[0,5],[5,4],[1,0],[4,1],[2,3],[4,2],[3,1],[5,2],[3,0]], "prism graph incidence changed")
    require([e["routing"] for e in definition["edges"]] == [["1","0","0","0"],["0","1","0","0"],["0","0","1","0"],["0","0","0","1"],["1","0","0","-1"],["0","1","0","-1"],["0","0","1","-1"],["1","-1","0","0"],["1","0","-1","0"]], "prism momentum routing changed")
    require(definition["targets"][1] == {"powers": [1,2,1,1,1,1,1,1,1], "numerator": "g1_2^2+g1_3*g2_4"} and definition["numerator_convention"] == "shifted_euclidean", "fixed original prism reference target changed")
    require(definition["chemical_potentials"] == ["1"] and all(e["mass_squared"] == "0" for e in definition["edges"]), "prism physical parameters changed")
    files.add(fixture)
    require((directory / "derivation-snapshot.md").is_file(), "saved independent derivation is missing")
    files.add(directory / "derivation-snapshot.md")
    expected_target = {**definition["targets"][1], "numerator_convention": definition["numerator_convention"]}
    proof_root = directory.parent / "prism-virtual-exact-reference"
    proof_hashes_path = proof_root / "artifact-hashes.json"
    proof_hashes = read(proof_hashes_path)["files"]
    files.add(proof_hashes_path)
    for name in ["status.json", "native-cargo-replay/input.json", "native-cargo-replay/result.json", "native-cargo-replay/native-reduction.json"]:
        path = proof_root / name
        require(sha(path) == proof_hashes[name], "native exact virtual proof artifact changed: " + name)
        files.add(path)
    proof_status = read(proof_root / "status.json")
    proof = read(proof_root / "native-cargo-replay/result.json")
    require(proof_status["candidate_fixture_sha256"] == sha(ROOT / proof_status["candidate_fixture"]), "virtual identity candidate fixture differs from exact proof")
    require(proof_status["runtime_required_external_solvers"] == [], "reference virtual identity verification cannot require an external solver")
    require([x["label"] for x in proof["relations"]] == ["C0A", "C1A", "C0B", "C1B", "CC"] and all(x["exact_zero"] and not x["residual_difference"] for x in proof["relations"]), "virtual coefficient identities lack complete exact source replay")
    require(proof["native_predictions_read"] == proof["oracle_records_read"] == 0, "virtual proof independence mismatch")
    bridge_root = directory.parent / "prism-support-wide-proof/two-virtual-thermal-bridge"
    bridge_status_path = bridge_root / "status.json"
    bridge_status = read(bridge_status_path)
    files.add(bridge_status_path)
    require(bridge_status["native_production_admission_changed"] is False, "reference proof cannot imply production admission")
    for name, info in bridge_status["files"].items():
        path = bridge_root / name
        require(sha(path) == info["sha256"], "thermal UV bridge proof artifact changed")
        files.add(path)
    bridge_audit = read(bridge_root / "two-virtual-uv-chart-audit.json")
    require(bridge_audit["status"] == "pass" and sum(r["supports"] for r in bridge_audit["runs"]) == 192, "thermal UV support audit incomplete")
    for run in bridge_audit["runs"]:
        path = bridge_root.parent / "native-full-kinematic" / Path(run["file"]).name
        require(sha(path) == run["sha256"], "exact routed support corpus changed")
        files.add(path)
    replay_path = directory / "evidence-replay.json"
    replay = read(replay_path)
    require(replay["status"] == "passed" and replay["sample_count"] == 60 and replay["exact_reports_replayed"] == 120, "exact regulator/tail evidence replay incomplete")
    require(replay["script_sha256"] == sha(directory / "replay_evidence.py"), "evidence replay source changed")
    files.update([replay_path, directory / "replay_evidence.py"])
    replay_rows = {(row["profile"], row["sample"]): row for row in replay["records"]}
    require(len(replay_rows) == 60, "duplicate or missing exact evidence replay rows")
    for name, digits, denominator, step, cutoff in PROFILES:
        location = directory / name
        metadata_path = location / "provenance.json"
        if not metadata_path.exists() or read(metadata_path)["status"] != "sample grid complete; refinements and Laurent acceptance separate":
            pending.append(name)
            continue
        metadata = read(metadata_path)
        files.add(metadata_path)
        config = metadata["configuration"]
        require((config["digits"], config["denominator"], config["step"], int(config["cutoff"]), config["count"]) == (digits, denominator, step, cutoff, 12), "reference configuration mismatch")
        require(metadata["oracle_answers_read"] == metadata["native_amf_predictions_read"] == 0, "reference independence mismatch")
        if reference_sources is None:
            reference_sources = metadata["source_sha256"]
            executable = metadata["executable_sha256"]
            for path, digest in reference_sources.items():
                path = ROOT / path
                require(sha(path) == digest, "reference generator source changed: " + str(path))
                files.add(path)
        require(metadata["source_sha256"] == reference_sources and metadata["executable_sha256"] == executable, "reference implementation changed across profiles")
        require(len(metadata["completed_samples"]) == 12, "incomplete grid metadata")
        samples = []
        for i in range(1, 13):
            stem = f"sample-{i:02}"
            path = location / (stem + ".json")
            row = read(path)
            files.add(path)
            require(row["digits"] == digits and row["step"] == step and row["truncation"] == cutoff, "sample configuration mismatch")
            require(Fraction(row["epsilon"]) == Fraction(i, denominator), "sample epsilon mismatch")
            require(row["target"] == expected_target, "reference target changed")
            require(row["two_cut_slots"] == [[1, 5], [1, 7], [5, 7]] and row["three_cut_slots"] == [1, 5, 7] and row["zero_cut_sectors"] == [[], [1], [5], [7]], "complete cut sum changed")
            require(len(row["rho_residue_checks"]) == 42 and all(c["passed"] for c in row["rho_residue_checks"]), "analytic regulator poles did not cancel")
            require(metadata["completed_samples"][i-1]["sha256"] == sha(path), "sample changed after generation")
            checks["analytic_regulator_residue_cancellations"] += 42
            regulator_path, tail_path = location/(stem+"-regulator.json"), location/(stem+"-tails.json")
            regulator, tail = read(regulator_path), read(tail_path)
            require(Fraction(regulator["epsilon"]) == Fraction(row["epsilon"]) == Fraction(tail["epsilon"]), "regulator/tail evidence epsilon differs from numerical sample")
            replay_row = replay_rows[(name, stem)]
            require(Fraction(replay_row["epsilon"]) == Fraction(row["epsilon"]), "replayed contour epsilon changed")
            for kind, evidence_path in [("regulator", regulator_path), ("tails", tail_path)]:
                checked = replay_row["checks"][kind]
                require(checked["replayed_report_equal"] and checked["saved_sha256"] == sha(evidence_path), "exact evidence differs from independently replayed report")
                require(checked["script_sha256"] == sha(directory / checked["script"]), "exact evidence analysis source changed")
            require(not regulator["problems"] and Fraction(regulator["minimum_nonconstant_gamma_pole_distance_at_rho_zero"]) > 0, "contour pole check failed")
            require(not tail["bad"] and Fraction(tail["minimum_gap"]) >= 2, "uniform tail damping failed")
            checks["sample_contour_pole_checks"] += 1
            checks["uniform_stirling_tail_checks"] += 1
            files.update([regulator_path, tail_path])
            archive_path = location/(stem+"-archives.json")
            archives = read(archive_path)
            files.add(archive_path)
            provenance_path = location/(stem+"-resources.provenance.json")
            provenance = read(provenance_path)
            files.add(provenance_path)
            require(provenance["argument_files"][config["executable"]]["sha256"] == executable, "sample executable provenance mismatch")
            for archive in archives:
                path = location/archive["file"]
                require(sha(path) == archive["sha256_gzip"], "compressed contour changed")
                require(hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest() == archive["sha256_uncompressed"], "lossless contour archive failed")
                files.add(path)
                checks["archived_contours_replayed"] += 1
                if archive["file"] == stem + "-contours.json.gz":
                    numerical_input = [info for key, info in provenance["argument_files"].items() if Path(key).name == stem + "-contours.json"]
                    require(len(numerical_input) == 1 and numerical_input[0]["sha256"] == archive["sha256_uncompressed"], "archived contour differs from actual numerical command input")
                    require(replay_row["contour_sha256"] == archive["sha256_uncompressed"] and replay_row["contour_archive_sha256"] == archive["sha256_gzip"], "pole/tail evidence was replayed against a different contour")
            samples.append(row)
        profiles.append({"name": name, "digits": digits, "epsilon_denominator": denominator, "integration_step": step, "integration_cutoff": cutoff, "fits": {field: interpolate(samples, field) for field in FIELDS}, "ten_node_interpolation_diagnostic": {field: interpolate(samples, field, 10) for field in FIELDS}})
    require(not pending or args.allow_incomplete, "missing completed independent profiles: " + ", ".join(pending))
    refinements = []
    for previous, current in zip(profiles, profiles[1:]):
        for field in FIELDS:
            for order in ORDERS:
                refinements.append({"from_profile": previous["name"], "to_profile": current["name"], "normalization_field": field, "order": order, **compare(previous["fits"][field][str(order)], current["fits"][field][str(order)])})
    raw_comparisons = [x for x in refinements if x["normalization_field"] == "full_raw_value"]
    order_checks = []
    imaginary_checks = []
    if profiles:
        finest_profile = profiles[-1]
        for order in ORDERS:
            order_checks.append({"order": order, **compare(finest_profile["fits"]["full_raw_value"][str(order)], finest_profile["ten_node_interpolation_diagnostic"]["full_raw_value"][str(order)])})
            imaginary = Decimal(finest_profile["fits"]["full_raw_value"][str(order)]["im"])
            imaginary_checks.append({"order": order, "absolute_imaginary_part": str(abs(imaginary)), "passed": abs(imaginary) <= Decimal("1e-25")})
    final_uncertainty = {}
    for order in ORDERS:
        changes = [Decimal(r["absolute_difference"]) for r in refinements if r["normalization_field"] == "full_msbar_reference_normalization" and r["order"] == order and r["from_profile"] != "laurent-baseline"]
        if profiles:
            fine = profiles[-1]
            diagnostic = compare(fine["fits"]["full_msbar_reference_normalization"][str(order)], fine["ten_node_interpolation_diagnostic"]["full_msbar_reference_normalization"][str(order)])
            changes.append(Decimal(diagnostic["absolute_difference"]))
        final_uncertainty[str(order)] = str(max(changes, default=Decimal(0)))
    final_raw_uncertainty = {}
    for order in ORDERS:
        changes = [Decimal(r["absolute_difference"]) for r in refinements if r["normalization_field"] == "full_raw_value" and r["order"] == order and r["from_profile"] != "laurent-baseline"]
        if profiles:
            fine = profiles[-1]
            diagnostic = compare(fine["fits"]["full_raw_value"][str(order)], fine["ten_node_interpolation_diagnostic"]["full_raw_value"][str(order)])
            changes.append(Decimal(diagnostic["absolute_difference"]))
        final_raw_uncertainty[str(order)] = str(max(changes, default=Decimal(0)))
    reference_uncertainty_passed = not pending and all(Decimal(v) <= Decimal("1e-12") for v in final_uncertainty.values())
    passed = not pending and len(raw_comparisons) == 20 and all(x["passed"] for x in raw_comparisons + order_checks + imaginary_checks) and reference_uncertainty_passed
    report = {
        "schema_version": 1,
        "status": "pending independent reference refinements" if pending else "passed" if passed else "failed",
        "scope": "complete supplemental raised original polynomial numerator prism reference; other prism target and native AMF comparisons remain separate",
        "target": expected_target,
        "working_pole_floor": -6,
        "requested_orders": [-4, 0],
        "normalization": "primary acceptance uses raw Euclidean coefficients; MSbar-normalized coefficients retained separately",
        "relative_tolerance": "1e-12", "small_magnitude_threshold": "1e-20", "small_absolute_tolerance": "1e-25",
        "uncertainty": "Empirical envelope of independently varied epsilon spacing, arithmetic precision, and quadrature step/cutoff; not a rigorous error bound. Small coefficients are retained as fitted values, not forced to zero.",
        "pending_profiles": pending,
        "target_absolute_reference_uncertainty_in_msbar_normalization": "1e-12",
        "final_empirical_absolute_uncertainty_in_msbar_normalization": final_uncertainty,
        "final_empirical_absolute_uncertainty_raw": final_raw_uncertainty,
        "reference_uncertainty_target_passed": reference_uncertainty_passed,
        "uncertainty_definition": "Maximum of the nearby/fine epsilon-grid change, fine precision change, fine quadrature change, and ten/twelve-node change. The initial coarse-grid change is retained separately as a convergence diagnostic, not used as an asserted bound on the final approximation.",
        "checks": checks,
        "exact_regulator_tail_reports_replayed": replay["exact_reports_replayed"],
        "virtual_identity_exact_checks": 5,
        "virtual_identity_scope": "Generic rational coefficient proof with retained nonzero conditions; its historical incomplete thermal/whole-reference flags are not reinterpreted as a complete reference theorem.",
        "executable_sha256": executable,
        "profiles": profiles,
        "refinements": refinements,
        "finest_interpolation_order_checks": order_checks,
        "finest_imaginary_zero_checks": imaginary_checks,
        "native_amf_predictions_compared": 0,
        "supplied_oracle_records_compared": 0,
        "arithmetic": "Decimal 130 digits, exact rational interpolation weights",
        "inputs_sha256": {relative_path(p): sha(p) for p in sorted(files)},
        "script_sha256": sha(Path(__file__)),
    }
    if profiles:
        finest = profiles[-1]["fits"]
        report["current_finest_coefficients"] = {field: {str(k): finest[field][str(k)] for k in ORDERS} for field in FIELDS}
        report["empirical_absolute_envelope"] = {field: {str(k): str(max((Decimal(r["absolute_difference"]) for r in refinements if r["normalization_field"] == field and r["order"] == k), default=Decimal(0))) for k in ORDERS} for field in FIELDS}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    if args.reference_output is not None:
        require(passed, "cannot publish a reference before all refinement gates pass")
        destination = args.reference_output.resolve()
        require(destination not in files and destination != output, "reference output cannot overwrite an input or validation report")
        values = report["current_finest_coefficients"]
        reference = {
            "schema_version": 1,
            "status": "independent_laurent_reference_empirically_validated",
            "definition": relative_path(fixture), "definition_sha256": sha(fixture),
            "target_index": 1, "target": expected_target,
            "dimension": "4-2*epsilon", "chemical_potential": "1",
            "complete_cut_subsets": [[], [1], [5], [7], [1,5], [1,7], [5,7], [1,5,7]],
            "normalization": "raw Euclidean amplitude",
            "laurent_orders": list(ORDERS),
            "coefficients": values["full_raw_value"],
            "empirical_absolute_uncertainty": final_raw_uncertainty,
            "comparison_normalization": "raw times (4*pi)^8*(exp(EulerGamma)/pi)^(4*epsilon)",
            "comparison_normalized_coefficients": values["full_msbar_reference_normalization"],
            "comparison_normalized_empirical_absolute_uncertainty": final_uncertainty,
            "uncertainty_caution": report["uncertainty"],
            "coefficients_forced_to_zero": False,
            "validation_report": relative_path(output), "validation_sha256": sha(output),
            "source_and_inputs_sha256": report["inputs_sha256"],
            "native_amf_predictions_read": 0, "supplied_oracle_records_read": 0,
            "native_four_loop_acceptance": False,
            "required_external_solvers": [],
        }
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(reference, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "pending": pending, "raw_refinement_comparisons": len(raw_comparisons), "passed_raw_comparisons": sum(r["passed"] for r in raw_comparisons)}))
    if not passed and not args.allow_incomplete:
        raise SystemExit(1)


if __name__ == "__main__":
    with localcontext() as context:
        context.prec = 130
        main()
