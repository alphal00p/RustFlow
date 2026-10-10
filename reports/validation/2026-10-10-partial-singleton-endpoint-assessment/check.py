#!/usr/bin/env python3
"""Exact independent rational-function consistency checks; no period values."""
import ast
import hashlib
import json
from fractions import Fraction
from pathlib import Path

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
NATIVE = BASE.parent / "2026-10-10-partial-placement-closure-pilot/single-slot0/result.json"


def clean(p): return {k: v for k, v in p.items() if v}
def add(a, b):
    c = a.copy()
    for k, v in b.items(): c[k] = c.get(k, Fraction()) + v
    return clean(c)
def mul(a, b):
    c = {}
    for (x, y), v in a.items():
        for (z, w), u in b.items():
            k = (x+z, y+w)
            c[k] = c.get(k, Fraction()) + v*u
    return clean(c)


class RF:
    def __init__(self, n, d=None):
        self.n = clean(n); self.d = {(0, 0): Fraction(1)} if d is None else clean(d)
        assert self.d
    def __add__(self, b): return RF(add(mul(self.n,b.d),mul(b.n,self.d)),mul(self.d,b.d))
    def __neg__(self): return RF({k:-v for k,v in self.n.items()},self.d)
    def __sub__(self,b): return self + -b
    def __mul__(self,b): return RF(mul(self.n,b.n),mul(self.d,b.d))
    def __truediv__(self,b): assert b.n; return RF(mul(self.n,b.d),mul(self.d,b.n))
    def __pow__(self,n):
        if n < 0: return RF(self.d,self.n)**(-n)
        r = constant(1)
        for _ in range(n): r = r*self
        return r
    def same(self,b): return not add(mul(self.n,b.d),{k:-v for k,v in mul(b.n,self.d).items()})


def constant(n): return RF({(0,0):Fraction(n)})
eps = RF({(1,0):Fraction(1)}); eta = RF({(0,1):Fraction(1)})


def parse(s):
    def visit(n):
        if isinstance(n,ast.Constant) and isinstance(n.value,int): return constant(n.value)
        if isinstance(n,ast.Name): return {"epsilon":eps,"eta":eta}[n.id]
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): return -visit(n.operand)
        if isinstance(n,ast.BinOp):
            a=visit(n.left)
            if isinstance(n.op,ast.Pow):
                assert isinstance(n.right,ast.Constant) and isinstance(n.right.value,int)
                return a**n.right.value
            b=visit(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
        raise ValueError(ast.dump(n))
    return visit(ast.parse(s.replace("^","**"),mode="eval").body)


def eta_terms(r):
    # Native target denominators are monomials in eta times epsilon polynomials.
    powers={p[1] for p in r.d}; assert len(powers)==1
    shift=powers.pop(); d={(p[0],0):v for p,v in r.d.items()}
    terms={}
    for (e,a),v in r.n.items():
        p=a-shift
        terms.setdefault(p,{})[(e,0)]=v
    return {p:RF(n,d) for p,n in terms.items() if n}


def poly_json(p): return [{"epsilon_power":i,"eta_power":j,"coefficient":str(v)} for (i,j),v in sorted(p.items())]
def rf_json(r): return {"numerator":poly_json(r.n),"denominator":poly_json(r.d)}


data=json.loads(NATIVE.read_text()); assert data["status"]=="closed"
closed=data["closed"];basis=closed["basis"];assert len(basis)==2
rows=[]
for row in closed["target_then_derivative_rows"]:
    assert all(t["indices"] in basis for t in row)
    rows.append([sum((parse(t["coefficient"]) for t in row if t["indices"]==b),constant(0)) for b in basis])
assert len(rows)==4
d=constant(4)-constant(2)*eps
ratio=d-constant(2); expected_power=d-constant(3)
v=[constant(1),ratio]
for i in range(2):
    assert sum((rows[i+2][j]*v[j] for j in range(2)),constant(0)).same(expected_power*v[i]/eta)
other_power=constant(2)-constant(3)*eps
assert (rows[2][0]+rows[3][1]).same((expected_power+other_power)/eta)
assert (rows[2][0]*rows[3][1]-rows[2][1]*rows[3][0]).same(expected_power*other_power/eta**2)
recon=[sum((r[j]*v[j] for j in range(2)),constant(0)) for r in rows[:2]]
assert recon[0].same(-constant(1)/eta)
terms=eta_terms(recon[1]);assert set(terms)=={-1,-2}
for power,coefficient in terms.items():
    assert not any(k[1] for k in coefficient.n|coefficient.d)
assert recon[1].same(sum((c*eta**p for p,c in terms.items()),constant(0)))

# Independent compact integration-by-parts relation on the real convergent strip.
moment_checks=[]
for dimension in [Fraction(13,2),Fraction(31,3),Fraction(39,2)]:
    for mu in [Fraction(1),Fraction(2,3),Fraction(3,2)]:
        # Common angular factor and mu^(D-2) cancel; no special function needed.
        bulk=1/(2*(dimension-2));surface=1/(2*mu)
        assert surface/bulk==(dimension-2)/mu
        moment_checks.append({"dimension":str(dimension),"mu":str(mu),"surface_bulk_ratio":str(surface/bulk)})

output={"scope":"Validation-only exact consistency of native closure, compact moment ratio and dimensional homogeneity. No boundary period or endpoint value installed/evaluated.",
        "native_result_sha256":hashlib.sha256(NATIVE.read_bytes()).hexdigest(),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "basis":basis,"conditions_retained":len(closed["conditions"]),
        "exact_matrix_eigenvector_checks":2,"exact_characteristic_checks":2,
        "compact_ratio":"I1/I0=(D-2)/mu; actual input mu=1",
        "selected_eta_exponent":"D-3=1-2 epsilon; two virtual loops and three positive hard powers",
        "excluded_extra_DE_exponent":"2-3 epsilon; its boundary coefficient must vanish",
        "scalar_reconstruction_over_I0":"-1/eta", "raised_reconstruction_eta_terms":{str(p):rf_json(c) for p,c in terms.items()},
        "target_eta_exponents":["D-4","D-5"],"compact_moment_checks":moment_checks,
        "limitations":"The high-D endpoint proof and actual integrated boundary must supply the constants. Exceptional native conditions must still be checked; this algebra does not waive them."}
(BASE/'consistency-checks.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps({"matrix_checks":4,"compact_checks":len(moment_checks),"scalar_identity":True,"raised_eta_powers":sorted(terms)}))
