"""Check saved singleton theorem certificates separately from native AMF proof.

This validates report structure/bindings and rational witness inequalities. It
is not a second theorem prover; the sealed owner is bound by the coherent build.
No independent reference value is read here.
"""
import json
from fractions import Fraction
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[3]/"tools/finite_density"))
from compare_massive_reference import complex_decimal

CONSTRUCTION="certified_massless_singleton_physical_zero"
PROOF="original-eta-independent-singleton-endpoint-zero-v1"
ORIGIN="joint-high-D-massless-flow-origin-v1;polynomial-completions;lower-contact-zero-jets;energy-regulator-removed-at-fixed-positive-T-and-eta;real-Fermi-T0;endpoint-then-meromorphic-D"
ENDPOINT="positive-eta-gaussian-germ-and-rank-one-blocks-v1;fixed-T-regulator-order;UV-meromorphic-before-high-D-endpoint"
CLASSIFICATIONS={"uniform-null-cone-endpoint","vanishing-required-cut","joint-dimensional-lower-contact-zero","scaleless-unrestricted-virtual-polynomial"}

def require(test,message):
    if not test: raise ValueError(message)

def exact_zero(text):
    return text=="0" or all(x==0 for x in complex_decimal(text))

def validate_certificates(records,definition,cuts):
    expected=[list(c) for c in cuts if len(c)==1]
    require(isinstance(records,list) and [r.get("cut_slots") for r in records]==expected,"exact singleton certificate cut list required")
    physical=len(definition["edges"])
    require(definition["loops"]>=2 and all(Fraction(e["mass_squared"])==0 for e in definition["edges"]),"certificate input is not wholly massless with a virtual loop")
    for item in records:
        cut=item["cut_slots"][0];c=item["certificate"]
        require(c["construction"]==CONSTRUCTION and c["proof_version"]==PROOF,"wrong physical-zero certificate owner")
        require(c["scope"]=="original eta-independent physical singleton targets only","wrong zero theorem scope")
        require(c["native_closure"] is False and c["transport_performed"] is False,"physical zero must not claim native closure or transport")
        require(all(isinstance(c[k],str) and c[k] for k in ("input_identity","proof_homotopy_identity","normalization")),"missing exact physical-zero family/input binding")
        require(isinstance(c["family_signature"],list) and c["family_signature"] and all(isinstance(x,str) and x for x in c["family_signature"]),"missing exact serialized family signature")
        require(c["proof_homotopy_shifted_slots"]==[i for i in range(physical) if i!=cut],"zero proof must use all uncut physical rows")
        require(isinstance(c["nonzero_conditions"],list) and all(isinstance(x,str) and x for x in c["nonzero_conditions"]),"missing retained coefficient domain")
        audit=c["label_audit"]; structure=audit["endpoint_structure"]
        require(audit["origin_identity"]==ORIGIN and structure["proof_version"]==ENDPOINT,"wrong origin/endpoint theorem identity")
        require(structure["virtual_loops"]==definition["loops"]-1,"singleton virtual rank differs from input")
        require([row[0] for row in structure["denominator_rows"]]==c["proof_homotopy_shifted_slots"],"incomplete singleton physical-row audit")
        require(structure["block_assignments"]==[] and structure["block_representatives"]==[] and structure["block_a_lower_bounds"]==[],"rank-one block evidence cannot masquerade as singleton germ")
        require(len(structure["block_q_upper_bounds"])==1 and Fraction(structure["block_q_upper_bounds"][0])>=0,"missing singleton Gaussian bound")
        labels=audit["labels"]
        require(isinstance(labels,list) and labels,"original target label audit missing")
        keys=[tuple(row["indices"]) for row in labels]
        require(len(keys)==len(set(keys)),"duplicate original target audit")
        for row in labels:
            require(row["classification"] in CLASSIFICATIONS,"non-singleton endpoint label")
            require(len(row["indices"])>=physical+2 and all(type(n)is int for n in row["indices"]),"bad physical label")
            d=Fraction(row["high_dimension_witness"])
            require(d>0,"invalid high-D witness")
            for bound in row["inequalities"]:
                degree=bound["degree"]
                require(Fraction(degree["twice_dimension_coefficient"],2)*d+degree["constant"]>0,"non-strict dimensional witness")
    return records

def load_certificates(directory,definition,cuts):
    path=directory/"physical-zero-certificates.json"
    return validate_certificates(json.loads(path.read_text()),definition,cuts),path

def validate_sample(record,certificates,definition,cuts,source_options=None):
    require(record.get("physical_zero_certificates")==certificates,"fixed prediction certificate array differs from saved preparation")
    occupied=record["occupied_reports"]
    require([r["cut_slots"] for r in occupied]==[list(c)for c in cuts[1:]],"occupied reports must retain all physical cuts")
    by_cut={tuple(r["cut_slots"]):r["certificate"]for r in certificates}
    contributions={tuple(r["cut_slots"]):r["values"]for r in record["contributions"]}
    for item in occupied:
        cut=tuple(item["cut_slots"]);evaluation=item.get("report",item)
        if cut in by_cut:
            c=by_cut[cut]
            require(evaluation["construction"]==CONSTRUCTION,"singleton must use the explicit physical-zero owner")
            require(evaluation["basis_size"]==0 and evaluation["native_storage_capacity"]is None and evaluation["source_options"]is None and evaluation["shifted_slots"]==[],"physical zero must not fabricate flow/storage/source metadata")
            require(evaluation["massless_endpoint"]==c["label_audit"] and evaluation["nonzero_conditions"]==c["nonzero_conditions"],"evaluated theorem/domain differs from preparation certificate")
            require(PROOF in evaluation["contour_admission"],"physical-zero evaluation lacks owner identity")
            require(evaluation["physical_arity"]==len(c["label_audit"]["labels"][0]["indices"]),"physical-zero arity differs from its audited labels")
            require(len(evaluation["values"])==len(definition["targets"]) and all(exact_zero(x)for x in evaluation["values"]),"physical-zero result must be exact zero for every original target")
            require(len(contributions[cut])==len(definition["targets"]) and all(exact_zero(x)for x in contributions[cut]),"assembled singleton contribution must match certified exact zeros")
        else:
            require(evaluation["construction"]=="weighted_amf" and evaluation["massless_endpoint"]is not None,"multicut requires actual native AMF and bound endpoint evidence")
            if source_options is not None:require(evaluation["source_options"]==source_options,"multicut source policy differs from saved configuration")
    require(record.get("vacuum_zero_certificates"),"missing native vacuum zero certificates")
    require(record.get("empty_support")==[],"no cut may be replaced by an empty-support assertion")
