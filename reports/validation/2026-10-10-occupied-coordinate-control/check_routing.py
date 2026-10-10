#!/usr/bin/env python3
"""Independent exact polynomial check; no native prediction or oracle input."""
import ast
import hashlib
import json
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def add(a,b,sign=1):
    out=defaultdict(Q,a)
    for k,v in b.items(): out[k]+=sign*v
    return {k:v for k,v in out.items() if v}


def mul(a,b):
    out=defaultdict(Q)
    for k,v in a.items():
        for j,w in b.items(): out[tuple(sorted(k+j))]+=v*w
    return {k:v for k,v in out.items() if v}


def poly(text):
    def visit(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int): return {():Q(n.value)}
        if isinstance(n,ast.Name): return {(n.id,):Q(1)}
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): return {k:-v for k,v in visit(n.operand).items()}
        if isinstance(n,ast.BinOp):
            a,b=visit(n.left),visit(n.right)
            if isinstance(n.op,ast.Add): return add(a,b)
            if isinstance(n.op,ast.Sub): return add(a,b,-1)
            if isinstance(n.op,ast.Mult): return mul(a,b)
        raise ValueError(ast.dump(n))
    return visit(ast.parse(text,mode='eval').body)


def main():
    out=HERE/'routing-audit.json';assert not out.exists()
    original_path=ROOT/'examples/finite_density/massless_three_loop_chain.json'
    new_path=HERE/'rerouted-input.json'
    a,b=json.loads(original_path.read_text()),json.loads(new_path.read_text())
    r=[[Q(1),Q(0),Q(-1)],[Q(0),Q(1),Q(-1)],[Q(0),Q(0),Q(-1)]]
    permutation=[3,4,2,0,1]
    assert sum(r[0][j]*(r[1][(j+1)%3]*r[2][(j+2)%3]-r[1][(j+2)%3]*r[2][(j+1)%3]) for j in range(3))==-1
    assert [[sum(r[i][k]*r[k][j] for k in range(3)) for j in range(3)] for i in range(3)]==[[Q(i==j) for j in range(3)] for i in range(3)]
    for new,old in enumerate(permutation):
        x,y=a['edges'][old],b['edges'][new]
        assert [sum(Q(x['routing'][k])*r[k][j] for k in range(3)) for j in range(3)]==list(map(Q,y['routing']))
        for key in ['vertices','charges','mass_squared']: assert x[key]==y[key]
    assert [[sum(r[i][j]*b['loop_charges'][j][species] for j in range(3)) for species in range(len(a['chemical_potentials']))] for i in range(3)]==a['loop_charges']
    assert a['chemical_potentials']==b['chemical_potentials']
    expected=defaultdict(Q)
    for i in range(3):
        for j in range(3):
            expected[(f'g{min(i,j)+1}_{max(i,j)+1}',)]+=r[0][i]*r[1][j]
            expected[tuple(sorted([f'u{i+1}',f'u{j+1}']))]+=r[0][i]*r[1][j]
    expected={k:v for k,v in expected.items() if v}
    assert poly(b['targets'][1]['numerator'])==expected
    assert poly(a['targets'][0]['numerator'])==poly(b['targets'][0]['numerator'])
    for x,y in zip(a['targets'],b['targets']): assert [x['powers'][i] for i in permutation]==y['powers']
    assert b['targets'][1]['powers'][3]==2
    assert list(map(Q,b['edges'][0]['routing']))==[1,0,0]
    report={
      'scope':'Exact rational routing/denominator/numerator/chemical/measure identity only; no production or numerical claim',
      'old_momenta_in_new':[[str(v) for v in row] for row in r],
      'old_edge_for_new':permutation,'determinant':-1,'absolute_jacobian':1,
      'all_five_edges_and_denominator_powers_match':True,'incidence_mass_charge_per_edge_unchanged':True,
      'loop_chemical_shifts_match':True,'original_raised_numerator_polynomial_matches':True,
      'raised_old_slot0_is_new_slot3':True,'selected_new_cut0_is_old_cut3':True,
      'normalization':'Unimodular change: absolute determinant one; cut orientation and chemical energy unchanged; no extra Wick or cut phase',
      'comparison_limit':'Combined routing, physical-edge order and completion-coordinate diagnostic, not an isolated completion change',
      'inputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [original_path,new_path,HERE/'run.py',Path(__file__)]}}
    out.write_text(json.dumps(report,indent=2)+'\n');print('exact routing audit PASS')


if __name__=='__main__':main()
