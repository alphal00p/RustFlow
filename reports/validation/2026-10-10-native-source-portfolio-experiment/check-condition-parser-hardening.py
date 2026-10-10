#!/usr/bin/env python3
"""Reject nonconstant exponents and verify archived comparison conclusions unchanged."""
import hashlib
import importlib.util
import json
from pathlib import Path

base=Path(__file__).resolve().parent
output=base/'condition-parser-hardening-check.json'
assert not output.exists()
spec=importlib.util.spec_from_file_location('exact_condition_v2',base/'condition-polynomials-v2.py')
helper=importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
bad=['eta^(2+epsilon)','eta^(epsilon+2)','(eta+1)^(1+epsilon)',
     'eta^epsilon','eta^(-1)','eta^(1/2)','eta^17','eta^(2+epsilon^2)']
for expression in bad:
    try:helper.parse(expression)
    except AssertionError:pass
    else:raise AssertionError('unsupported exponent admitted: '+expression)
good=['eta^0','eta^2','eta^(4/2)','eta^(2+epsilon-epsilon)']
for expression in good:helper.parse(expression)
assert helper.parse('eta^(2+epsilon-epsilon)')==helper.parse('eta^2')
files=[base/'condition-parser-hardening-provenance.json',base/'condition-polynomials.py',
       base/'condition-polynomials-v2.py',base/'compare-portfolio.py',
       base/'compare-portfolio-hardened.py',Path(__file__)]
comparisons=[]
for label in ['original-19','strongest62']:
    old_path=base/(label+'-portfolio-comparison.json')
    new_path=base/(label+'-portfolio-comparison-hardened.json')
    old=json.loads(old_path.read_text());new=json.loads(new_path.read_text())
    assert {k:v for k,v in old.items()if k!='inputs_sha256'}=={k:v for k,v in new.items()if k!='inputs_sha256'}
    comparisons.append({'label':label,'all_saved_conclusions_rows_conditions_and_proof_hashes_identical':True,
        'points':new['points'],'condition_gate_pass_points':new['condition_gate_pass_points']})
    files.extend([old_path,new_path])
report={'status':'passed','scope':'Validation-only parser hardening. Archived conditions use supported literal powers; reruns give identical rows, parameter-coverage decisions, proof hashes and counts. Original files preserved.',
        'rejected_exponent_controls':bad,'accepted_exponent_controls':good,
        'comparison_reruns':comparisons,
        'artifact_sha256':{path.name:hashlib.sha256(path.read_bytes()).hexdigest()for path in files},
        'production_modified':False,'new_numerical_runs':0}
output.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'rejected':len(bad),'accepted':len(good),'identical_comparison_reruns':len(comparisons)}))
