#!/usr/bin/env python3
"""Exact finite-degree raw-C1 Ward tests, independent of native predictions."""
import hashlib
import json
from fractions import Fraction as Q
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / '2026-10-10-factorized-compact-ward-assessment/check.py'
spec = importlib.util.spec_from_file_location('moment_check', BASE)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
moment = module.moment


def main():
    out = HERE/'checks.json'; assert not out.exists()
    rows=[]
    for d in [7,9,11]:
        for k in range(7):
            for s in range(5):
                for c in [Q(-3,2),Q(2,3)]:
                    mu=Q(5,4)
                    m=c**k*moment(d,1,s,k,mu)
                    surface=c**k*moment(d,1,s+1,k+1,mu)
                    residual=(d-2+k)*m+(-1 if s==0 else s)*surface
                    assert residual==0
                    rows.append({'D':d,'energy_degree':k,'upper':s,'factor_coefficient':str(c),'residual':'0'})
    # A Gram-polynomial off-cone test function makes the higher-cut exclusion
    # explicit: F(q^2)=1+q^2, q^2 C2=C1. Its two homogeneous pieces cannot share
    # the D-4+k coefficient that works for a genuinely factorized C2 shell.
    d,k,mu=9,2,Q(5,4)
    m=moment(d,2,0,k,mu)+moment(d,1,0,k,mu)
    surface=moment(d,2,1,k+1,mu)+moment(d,1,1,k+1,mu)
    counter=(d-4+k)*m-surface
    assert counter==-2*moment(d,1,0,k,mu) and counter!=0
    # The old C1 coefficient is insufficient as soon as k>0.
    old=(d-2)*moment(d,1,0,k,mu)-moment(d,1,1,k+1,mu)
    assert old==-k*moment(d,1,0,k,mu) and old!=0
    report={'scope':'Exact raw-C1 distribution moment and guard-counterexample checks only; no source policy, closure or numerical prediction',
        'checks':len(rows),'samples':rows,'old_diagonal_counterexample':str(old),
        'C2_off_cone_gram_factor_counterexample':str(counter),
        'native_or_oracle_values_read':0,
        'inputs_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE,Path(__file__)]}}
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'checks':len(rows),'status':'passed','counterexamples_nonzero':True}))


if __name__=='__main__':main()
