"""Compare saved exploratory RustFlow transports with independent quadrature."""
import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

getcontext().prec = 100
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--predictions", type=Path, required=True)
parser.add_argument("--reference", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
predictions = json.loads(args.predictions.read_text())
reference = json.loads(args.reference.read_text())
# Bind the full exploratory input and every claimed refinement configuration.
assert predictions["schema"] == reference["schema"] == 1
assert predictions["status"] == "predictions saved before comparison"
assert predictions["reference_values_read"] is False
assert reference["status"] == "passed"
assert reference["scope"] == "validation-only exploratory one-loop flow and scalar mixed sunset asymptotic; no production acceptance"
assert reference["production_reference_values_loaded"] is False
assert reference["four_loop_acceptance"] is False
assert predictions["normalization"] == "angular volume omitted consistently"
assert reference["angular_normalization"] == "radial integrals omit Omega_(d-1)/(2*pi)^d; the virtual bubble retains (4*pi)^(-D/2)"
assert [Fraction(predictions[k]) for k in ("mass_squared", "chemical_potential", "radius_squared")] == [Fraction(1,4), Fraction(1), Fraction(3,4)]
dimensions = {Fraction(7,5), Fraction(3), Fraction(9,2)}
profiles = {(60,48,8,32), (60,64,8,48), (80,64,8,48), (80,64,16,48), (80,64,16,64)}
paths = {"positive-real", "upper-half-plane"}
expected_grid = {(d, *profile, path) for d in dimensions for profile in profiles for path in paths}
actual_grid = [(Fraction(row["spatial_dimension"]), *(row[k] for k in
    ("working_digits", "transport_order", "start_eta", "boundary_order", "path"))) for row in predictions["records"]]
assert len(actual_grid) == len(set(actual_grid)) == 30
assert set(actual_grid) == expected_grid

def finite_complex(value):
    re, im = Decimal(value["re"]), Decimal(value["im"])
    assert re.is_finite() and im.is_finite()
    return re, im

reference_grid = set()
for row in reference["one_loop"]:
    re, im = finite_complex(row["eta"])
    if (re, im) == (0, 0):
        eta_key = "zero"
    elif (re, im) == (16, 0):
        eta_key = "sixteen"
    elif (re, im) == (1, Decimal("0.5")):
        eta_key = "complex"
    else:
        assert im == 0 and abs(3*re-1) < Decimal("1e-55")
        eta_key = "one-third"
    assert finite_complex(row["mass_squared"]) == (Decimal("0.25"), 0)
    assert finite_complex(row["radius_squared"]) == (Decimal("0.75"), 0)
    key = (row["digits"], Fraction(row["spatial_dimension"]), eta_key)
    assert key not in reference_grid
    reference_grid.add(key)
assert reference_grid == {(digits,d,eta) for digits in (60,90) for d in dimensions
    for eta in ("zero", "sixteen", "complex", "one-third")}
mass_grid = set()
for row in reference["physical_mass_derivatives"]:
    key = (row["digits"], Fraction(row["spatial_dimension"]))
    assert key not in mass_grid
    mass_grid.add(key)
    assert row["original_numerator"] == "g+u^2; g=q^2 Minkowski, u=q^0"
    assert row["fixed_input_coefficient_eta_extension"] == "a*I + K, K=int E_eta/2"
    assert row["eta_only_raising_is_incorrect"] is True
assert mass_grid == {(digits,d) for digits in (60,90) for d in dimensions}
refs = {}
for record in reference["one_loop"]:
    if record["digits"] == 90 and Decimal(record["eta"]["re"]) == 0 and Decimal(record["eta"]["im"]) == 0:
        assert Fraction(record["spatial_dimension"]) not in refs
        refs[Fraction(record["spatial_dimension"])] = {"I":record["I"], "S":record["S"], "eta_derivative":record["I_prime"]}
for record in reference["physical_mass_derivatives"]:
    if record["digits"] == 90:
        refs[Fraction(record["spatial_dimension"])].update({key:record[key] for key in
            ("physical_scalar_mass_derivative", "physical_medium_mass_derivative")})
comparisons=[]
refinements=[]
for record in predictions["records"]:
    dimension = Fraction(record["spatial_dimension"])
    for key,expected in refs[dimension].items():
        finite_complex(record[key])
        finite_complex(expected)
        assert Decimal(expected["re"]) != 0
        relative=abs(Decimal(record[key]["re"])-Decimal(expected["re"]))/abs(Decimal(expected["re"]))
        imaginary=abs(Decimal(record[key]["im"])-Decimal(expected["im"]))
        assert relative < Decimal("1e-28"), (record,key,relative)
        assert imaginary < Decimal("1e-28"), (record,key,imaginary)
        comparisons.append({"dimension":str(dimension),"quantity":key,"configuration":[record[k] for k in
            ("working_digits","transport_order","start_eta","boundary_order","path")],
            "relative_error":str(relative),"imaginary_absolute_error":str(imaginary)})
    finest = next(r for r in predictions["records"] if Fraction(r["spatial_dimension"]) == dimension
        and r["working_digits"] == 80 and r["start_eta"] == 16 and r["boundary_order"] == 64 and r["path"] == "positive-real")
    if record is not finest:
        for key in refs[dimension]:
            relative=abs(Decimal(record[key]["re"])-Decimal(finest[key]["re"]))/abs(Decimal(finest[key]["re"]))
            assert relative < Decimal("1e-28")
            refinements.append({"dimension":str(dimension),"quantity":key,"relative_error":str(relative)})
result={"schema":1,"status":"passed","scope":"exploratory radial ODE and fixed original numerator mass derivatives, not production finite-density acceptance",
    "validated_contract":{"exact_unique_transport_grid":30,"dimensions":["7/5","3","9/2"],"profiles":[list(x) for x in sorted(profiles)],"paths":sorted(paths),"unique_reference_quadratures":24,"unique_mass_derivative_checks":6,"supported_schema":1},
    "relative_tolerance":"1e-28","imaginary_absolute_tolerance":"1e-28",
    "reference_comparisons":len(comparisons),"refinement_comparisons":len(refinements),
    "maximum_relative_reference_error":str(max(Decimal(x["relative_error"]) for x in comparisons)),
    "maximum_relative_refinement_error":str(max(Decimal(x["relative_error"]) for x in refinements)),
    "inputs_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (args.predictions,args.reference,Path(__file__))},"comparisons":comparisons,"refinements":refinements}
args.output.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:v for k,v in result.items() if k not in ("comparisons","refinements","inputs_sha256")}))
