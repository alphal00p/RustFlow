#!/usr/bin/env python3
"""Compare completed native three-loop predictions to separately generated references.

Reads saved artifacts only. Requires four sample or five Laurent profiles, with
independent precision/order/start/grid changes. No production solver is imported.
"""
import argparse
import json
import re
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/"tools/finite_density"))
from certificate_validation import load_certificates, validate_sample
from compare_massive_reference import (
    ROOT, ZERO, CONFIGURATIONS, complex_decimal, difference, digest, read,
)

TARGETS = ("scalar", "raised_original_numerator")
CUTS = ((), (0,), (3,), (0, 3))
SECTORS = ("vacuum", "single_cut_0", "single_cut_3", "double_cut", "total")
SETTINGS = ("digits", "series_order", "occupied_start_scale", "epsilon_grid_denominator")


def require(value, message):
    if not value:
        raise ValueError(message)


def label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def value(text):
    return (ZERO, ZERO) if text == "0" else complex_decimal(text)


def reference_adapter(path):
    raw = read(path)
    validation_path = path.with_name("validation.json")
    source_path = path.with_name("source-sha256.json")
    validation, sources = read(validation_path), read(source_path)
    require(validation["reference_sha256"] == digest(path), "independent reference changed after validation")
    require(validation["status"].startswith("passed independent reference"), "reference validation incomplete")
    require(raw["oracle_answers_read"] == 0 and raw["native_amf_predictions_read"] == 0, "reference independence mismatch")
    require(raw["normalization"] == "raw Euclidean", "reference normalization mismatch")
    require(raw["production_shortcuts_added"] == 0, "reference cannot be a production shortcut")
    inputs = [path, validation_path, source_path]
    for name, sha in sources.items():
        source = ROOT / name
        require(digest(source) == sha, "independent reference source artifact changed: " + name)
        inputs.append(source)
    samples = []
    for row in raw["samples"]:
        if row["digits"] != 80:
            continue
        require(tuple(tuple(c["cuts"]) for c in row["contributions"]) == CUTS, "independent reference cut list mismatch")
        samples.append({"dimension": str(Fraction(4)-2*Fraction(row["epsilon"])), "chemical_potential": "1", "total": row["values"], **{sector: c["values"] for sector,c in zip(SECTORS[:-1],row["contributions"])}})
    expansions = [r for r in raw["laurent"] if r["digits"] == 80]
    require(len(expansions) == 1, "independent Laurent reference is not unique")
    return {"status":"independent_reference_generated", "normalization":"raw_euclidean_dDp_over_2pi_power_D_per_loop", "oracle_records_read":0,"native_predictions_read":0,"source":{"input":raw["definition"],"input_sha256":sources[raw["definition"]]},"samples":samples,"laurent":{"orders":expansions[0]["orders"],"chemical_potential":"1","total":expansions[0]["values"]},"uncertainty":validation["accuracy"]}, inputs


def validate_native_provenance(args, directories, definition_path):
    # The runtime harness prepares all sectors through Run::full in both modes.
    # Laurent records currently omit per-sector construction reports, so bind
    # them to a completed fixed-D gate using the same input, source options,
    # closure options, exact executable and native compile identities.
    provenance_path = args.run_provenance.resolve()
    snapshot_path = args.source_snapshot.resolve()
    provenance_paths = [provenance_path] + [p.resolve() for p in args.additional_run_provenance]
    require(len(provenance_paths) == len(directories), "one runtime provenance per prediction directory required")
    snapshot = read(snapshot_path)
    executable = None
    for captured in provenance_paths:
        provenance = read(captured)
        executables = [v for k,v in provenance["argument_files"].items() if Path(k).name.startswith("finite_density_runtime_flow-")]
        require(len(executables) == 1, "no unique captured runtime native executable")
        sha = executables[0]["sha256"]
        require(executable is None or executable == sha, "native executable changed between resumed runs")
        executable = sha
    inputs = provenance_paths + [snapshot_path]
    checked_sources = {}
    for name in ["Cargo.lock","src/finite_density/assembly.rs","src/finite_density/flow.rs","src/finite_density/guarded.rs","src/finite_density/singleton_zero.rs","src/finite_density/massless_endpoint.rs","tests/finite_density_runtime_flow.rs",str(definition_path.relative_to(ROOT))]:
        require(name in snapshot["files"], "native source snapshot missing " + name)
        require(digest(ROOT/name) == snapshot["files"][name], "native source artifact changed: " + name)
        checked_sources[name] = snapshot["files"][name]
    baseline = None
    compile_ids = set()
    metadata_files = []
    certificate_arrays = []
    for directory in directories:
        certs, certpath = load_certificates(directory,read(definition_path),CUTS)
        certificate_arrays.append(certs); inputs.append(certpath)
        config_path = directory/"configuration.json"
        config = read(config_path);inputs.append(config_path)
        require(config["independent_reference_comparisons"] == 0, "native configuration contains reference comparisons")
        require(config["input_path"] == str(definition_path.relative_to(ROOT)), "native configuration uses different input")
        require(config["profiles"] == [list(x) for x in CONFIGURATIONS], "full precision/order/start profiles required")
        # Rust Debug represents checkpoint paths as escaped quoted strings.
        closure = re.sub(r'checkpoints: Some\((?:\\"|\").*?(?:\\"|\")\)', 'checkpoints: None', config["native_closure"])
        normalized = {"source_options":config["source_options"],"native_closure":closure,"guard_digits":config["guard_digits"]}
        require(baseline is None or normalized == baseline, "native source/closure policy changed between prediction directories")
        baseline = normalized
        for cuts in (CUTS[-1],):
            sub = directory/"native-closure"/("cut-"+"-".join(map(str,cuts)))
            files = sorted(sub.glob("round-*-closed.json"))
            require(files, "missing successful native closure metadata for cut " + str(cuts))
            for path in files:
                row=read(path);require(row["status"]=="closed","unclosed native program")
                fields=row["measure_id"].split(":",3)
                require(fields[0]=="rustflow-weighted-sources-v1" and all(re.fullmatch(r"[0-9a-f]{64}",x) for x in fields[1:3]),"missing native dependency compile identities")
                compile_ids.add(tuple(fields[1:3]));inputs.append(path);metadata_files.append({"path":label(path),"sha256":digest(path),"cut_slots":cuts,"physical_arity":row["physical_arity"],"storage_capacity":row["storage_capacity"]})
    require(len(compile_ids)==1,"native dependency compile identities changed")
    require(all(c==certificate_arrays[0]for c in certificate_arrays),"physical singleton certificates changed between resumed runs")
    report={"configuration":baseline,"executable_sha256":executable,"dependency_compile_ids":list(next(iter(compile_ids))),"source_snapshot_sha256":digest(snapshot_path),"checked_source_files":checked_sources,"closed_metadata":metadata_files,"physical_zero_certificates":certificate_arrays[0]}
    if args.mode=="laurent":
        require(args.sample_comparison is not None,"Laurent gate requires matching successful fixed-D construction evidence")
        sample_path=args.sample_comparison.resolve();sample=read(sample_path);inputs.append(sample_path)
        require(sample["status"]=="passed" and sample["mode"]=="sample","fixed-D construction evidence did not pass")
        old=sample["native_provenance"]
        for key in ["configuration","executable_sha256","dependency_compile_ids","source_snapshot_sha256","physical_zero_certificates"]:
            require(old[key]==report[key],"Laurent provenance differs from successful fixed-D run: "+key)
        for name,sha in sample["inputs_sha256"].items():
            path=ROOT/name if not Path(name).is_absolute() else Path(name)
            require(digest(path)==sha,"fixed-D supporting artifact changed: "+name)
        report["fixed_dimension_comparison"]={"path":label(sample_path),"sha256":digest(sample_path)}
        report["construction_limit"]="Laurent harness omits direct per-sector reports; multicut weighted-AMF construction and separate singleton physical-zero certificates are bound to the successful fixed-D gate and the identical Run::full source/executable and exact native compile identities. No per-sector Laurent coefficient claim."
    return report,inputs

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("sample", "laurent"))
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--additional-predictions", type=Path, action="append", default=[])
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--run-provenance", type=Path, required=True)
    parser.add_argument("--additional-run-provenance", type=Path, action="append", default=[])
    parser.add_argument("--source-snapshot", type=Path, required=True)
    parser.add_argument("--sample-comparison", type=Path, help="required for Laurent: successful matching full fixed-D comparison")
    args = parser.parse_args()
    laurent = args.mode == "laurent"
    dirs = [args.predictions.resolve()] + [p.resolve() for p in args.additional_predictions]
    require(len(set(dirs)) == len(dirs), "duplicate prediction directory")
    refpath = args.reference.resolve()
    reference, reference_inputs = reference_adapter(refpath)
    require(reference["status"] == "independent_reference_generated", "reference generation incomplete")
    require(reference["normalization"] == "raw_euclidean_dDp_over_2pi_power_D_per_loop", "reference normalization mismatch")
    require(reference["oracle_records_read"] == 0 and reference["native_predictions_read"] == 0, "reference must be independent")
    definition_path = ROOT / reference["source"]["input"]
    require(digest(definition_path) == reference["source"]["input_sha256"], "reference input digest changed")
    definition = read(definition_path)
    require(definition["name"] == "massless_three_loop_chain", "unexpected reference definition")
    input_paths = [d / "input.json" for d in dirs]
    require(all(read(p) == definition for p in input_paths), "prediction definition differs from independent reference")
    configs = tuple((*c, 1000) for c in CONFIGURATIONS) + ((28, 80, 12, 2000),) if laurent else tuple(CONFIGURATIONS)
    paths = []
    duplicate_paths = []
    for config in configs:
        digits, order, start = config[:3]
        suffix = f"-grid-{config[3]}" if laurent else ""
        name = f"prediction-{digits}-{order}-{start}{suffix}.json"
        candidates = [d / name for d in dirs if (d / name).is_file()]
        require(candidates, f"missing completed native profile {config}")
        require(all(digest(p) == digest(candidates[0]) for p in candidates), "conflicting duplicated profile")
        paths.append(candidates[0]); duplicate_paths.extend(candidates[1:])
    inputs = reference_inputs + [definition_path] + input_paths + paths + duplicate_paths
    native_provenance, native_inputs = validate_native_provenance(args, dirs, definition_path)
    inputs.extend(native_inputs)
    output = args.output.resolve()
    require(output not in inputs, "output cannot replace an input artifact")
    certificate_arrays=[]
    for directory in dirs:
        certificates, certificate_path=load_certificates(directory,definition,CUTS)
        certificate_arrays.append(certificates); inputs.append(certificate_path)
    require(all(c==certificate_arrays[0]for c in certificate_arrays),"physical-zero certificates changed across resumed runs")
    inputs.append(Path(__file__).with_name("certificate_validation.py"))
    previous = None
    guard = None
    dimension = None
    comparisons = []
    refinements = []
    assembly = []
    for profile_index, (config, path) in enumerate(zip(configs, paths)):
        record = read(path)
        require(record["full_amplitude"] is True, "complete amplitude required")
        require(record["normalization"] == "unscaled Euclidean amplitude", "native normalization mismatch")
        require(record["independent_reference_comparisons"] == 0, "native outputs must precede comparison")
        require(tuple(record[k] for k in SETTINGS[:3]) == config[:3], "native configuration mismatch")
        require(type(record["guard_digits"]) is int and record["guard_digits"] >= 0, "invalid guard precision")
        guard = record["guard_digits"] if guard is None else guard
        require(record["guard_digits"] == guard, "guard precision changed between profiles")
        require(record["guard_digits"] == native_provenance["configuration"]["guard_digits"], "prediction guard precision differs from saved native configuration")
        values = {}
        expected = {}
        if laurent:
            require(record["epsilon_grid_denominator"] == config[3], "epsilon-grid mismatch")
            require(len(record["expansions"]) == len(TARGETS), "target count mismatch")
            require(reference["laurent"]["orders"] == [-3, -2, -1, 0], "reference coefficient orders mismatch")
            require(reference["laurent"]["chemical_potential"] == definition["chemical_potentials"][0], "reference chemical potential mismatch")
            for i, expansion in enumerate(record["expansions"]):
                require(type(expansion["verified_digits"]) is int and expansion["verified_digits"] >= config[0], "shared Laurent fitter did not verify requested digits")
                require(set(expansion["coefficients"]) == {"-3", "-2", "-1", "0"}, "requested coefficient set mismatch")
                for j, power in enumerate((-3, -2, -1, 0)):
                    key = (TARGETS[i], "total", power)
                    values[key] = value(expansion["coefficients"][str(power)])
                    expected[key] = value(reference["laurent"]["total"][i][j])
        else:
            d = Fraction(4) - 2 * Fraction(record["epsilon"])
            dimension = d if dimension is None else dimension
            require(d == dimension, "dimension changed between profiles")
            candidates = [s for s in reference["samples"] if Fraction(s["dimension"]) == d and Fraction(s["chemical_potential"]) == Fraction(definition["chemical_potentials"][0])]
            require(len(candidates) == 1, "no unique independent reference for sample")
            sample = candidates[0]
            require(tuple(tuple(c["cut_slots"]) for c in record["contributions"]) == CUTS, "complete physical cut set mismatch")
            occupied = record["occupied_reports"]
            require(tuple(tuple(c["cut_slots"]) for c in occupied) == CUTS[1:], "occupied provenance does not cover all cuts")
            validate_sample(record, certificate_arrays[dirs.index(path.parent)], definition, CUTS,
                            read(path.parent/"configuration.json")["source_options"])
            require(record.get("vacuum_zero_certificates"), "missing native vacuum zero certificates")
            require(record.get("empty_support") == [], "no occupied sector may be replaced by an empty-support assertion in this gate")
            for sector, entries in zip(SECTORS, [c["values"] for c in record["contributions"]] + [record["values"]]):
                require(len(entries) == len(TARGETS), "target count mismatch")
                for i, target in enumerate(TARGETS):
                    key = (target, sector, None)
                    values[key] = value(entries[i]); expected[key] = value(sample[sector][i])
            for target in TARGETS:
                summed = tuple(sum((values[(target, sector, None)][part] for sector in SECTORS[:-1]), ZERO) for part in range(2))
                assembly.append({"target": target, "configuration": config, **difference(summed, values[(target, "total", None)])})
        for key, actual in values.items():
            target, sector, power = key
            comparisons.append({"target": target, "sector": sector, "power": power, "configuration": config, "prediction_path": label(path), "prediction_real": str(actual[0]), "prediction_imaginary": str(actual[1]), "reference_real": str(expected[key][0]), **difference(actual, expected[key]), "imaginary_zero_passed": abs(actual[1]) <= Decimal("1e-25")})
            if previous is not None:
                refinements.append({"target": target, "sector": sector, "power": power, "from_configuration": configs[profile_index-1], "to_configuration": config, "varied_setting": SETTINGS[profile_index-1], **difference(actual, previous[key])})
        previous = values
    passed = all(x["passed"] and x["imaginary_zero_passed"] for x in comparisons) and all(x["passed"] for x in refinements + assembly)
    report = {"schema_version": 1, "status": "passed" if passed else "failed", "mode": args.mode, "scope": "complete nonfactorized massless three-loop chain: separate certified physical singleton zeros plus native AMF double cut; no four-loop acceptance", "normalization": "unscaled Euclidean amplitude", "dimension": str(dimension) if dimension else "4-2*epsilon", "configurations": configs, "guard_digits": guard, "native_provenance": native_provenance, "relative_tolerance": "1e-12", "small_magnitude_threshold": "1e-20", "small_absolute_tolerance": "1e-25", "imaginary_zero_absolute_tolerance": "1e-25", "reference_comparison_count": len(comparisons), "independent_refinement_comparison_count": len(refinements), "supplied_oracle_numerical_records_compared": 0, "reference_uncertainty": reference["uncertainty"], "arithmetic": "Python Decimal, 100 digits, saved decimal inputs", "inputs_sha256": {label(p): digest(p) for p in inputs}, "script_sha256": digest(Path(__file__)), "comparisons": comparisons, "refinements": refinements, "assembly_checks": assembly}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("status", "mode", "reference_comparison_count", "independent_refinement_comparison_count")}))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
