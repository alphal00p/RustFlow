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
assert predictions["reference_values_read"] is False
assert reference["status"] == "passed"
assert len(predictions["records"]) == 30
refs = {}
for record in reference["one_loop"]:
    if record["digits"] == 90 and Decimal(record["eta"]["re"]) == 0 and Decimal(record["eta"]["im"]) == 0:
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
    "relative_tolerance":"1e-28","imaginary_absolute_tolerance":"1e-28",
    "reference_comparisons":len(comparisons),"refinement_comparisons":len(refinements),
    "maximum_relative_reference_error":str(max(Decimal(x["relative_error"]) for x in comparisons)),
    "maximum_relative_refinement_error":str(max(Decimal(x["relative_error"]) for x in refinements)),
    "inputs_sha256":{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        (args.predictions,args.reference,Path(__file__))},"comparisons":comparisons,"refinements":refinements}
args.output.write_text(json.dumps(result,indent=2)+"\n")
print(json.dumps({k:v for k,v in result.items() if k not in ("comparisons","refinements","inputs_sha256")}))
