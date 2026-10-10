#!/usr/bin/env python3
"""Source-descriptor and exact generator-algebra audit; not native rule replay."""
import ast
import hashlib
import json
import re
import struct
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from check import inverse, moment

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
META = ROOT / "reports/validation/2026-10-10-requested-ray-point-integration/requested-ray-point-three-loop-double/native-closure/round-000-provisional.json"


class Reader:
    def __init__(self, data): self.data, self.pos = data, 0
    def take(self, n):
        assert 0 <= n <= len(self.data) - self.pos
        result = self.data[self.pos:self.pos+n]; self.pos += n; return result
    def uint(self):
        x = self.take(1)[0]
        if x < 251: return x
        assert x in [251, 252, 253, 254]
        return int.from_bytes(self.take(2 ** (x-250)), "little")
    def sint(self):
        x = self.uint(); return x//2 if x % 2 == 0 else -(x+1)//2
    def option(self):
        x = self.take(1)[0]; assert x in [0, 1]
        return self.sint() if x else None
    def vec(self, f): return [f() for _ in range(self.uint())]
    def string(self): return self.take(self.uint()).decode()
    def domain(self): return self.vec(lambda: [self.option(), self.option()])
    def label(self): return self.vec(lambda: [self.take(1)[0], self.sint()])
    def term(self): return {"integral": self.label(), "coefficient_id": self.uint()}
    def source(self):
        return {"id": self.string(), "domain": self.domain(),
                "condition_ids": self.vec(self.uint), "terms": self.vec(self.term)}


def affine(text):
    """Read the actual saved affine factors, without evaluating arbitrary code."""
    text = re.sub(r"rustflow_occupied::q_(\d+)", r"q\1", text)
    def visit(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return [Q(0)] * 9 + [Q(node.value)]
        if isinstance(node, ast.Name) and re.fullmatch(r"q[0-8]", node.id):
            return [Q(i == int(node.id[1:])) for i in range(9)] + [Q(0)]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return [-x for x in visit(node.operand)]
        if isinstance(node, ast.BinOp):
            a, b = visit(node.left), visit(node.right)
            if isinstance(node.op, ast.Add): return [x+y for x,y in zip(a,b)]
            if isinstance(node.op, ast.Sub): return [x-y for x,y in zip(a,b)]
            if isinstance(node.op, ast.Mult):
                if any(a[:9]): a,b = b,a
                assert not any(a[:9]); return [a[9]*x for x in b]
        raise ValueError(ast.dump(node))
    return visit(ast.parse(text, mode="eval").body)


def main():
    out = HERE / "redundancy-checks.json"
    assert not out.exists(), "completed evidence is immutable"
    metadata = json.loads(META.read_text())
    raw = META.with_suffix(".bin").read_bytes()
    assert raw[:8] == b"RRPBIN\r\n"
    assert int.from_bytes(raw[8:12], "little") == 1
    pos = 20; sections = {}
    for _ in range(int.from_bytes(raw[16:20], "little")):
        tag, reserved, n = struct.unpack_from("<HHQ", raw, pos)
        assert reserved == 0; pos += 12
        sections[tag] = raw[pos:pos+n]; pos += n
    assert pos == len(raw)
    rd = Reader(sections[4])
    assert rd.string() == "rustred.guarded-source-program.v2"
    measure_id = rd.string(); assert measure_id == metadata["measure_id"]
    roles = rd.vec(lambda: rd.take(1)[0]); indices = rd.vec(rd.uint)
    sources = rd.vec(rd.source); zeros = rd.vec(rd.domain)
    # Deliberately no coefficient import, generated-rule decoding, or replay claim.
    payload, _ = json.JSONDecoder().raw_decode(measure_id.split(":", 3)[3])
    factor_text = payload["measure"].split("; factors=[", 1)[1].split("];", 1)[0]
    factors = [affine(x.strip()) for x in factor_text.split(",")]
    assert len(factors) == 13 and len(roles) == 16
    assert roles == [1,0,0,1,0,0,0,0,0,2,2,2,2,0,0,0]
    mat = [row[:9] for row in factors[:9]]; inv = inverse(mat)
    pairs = [(0,0),(0,1),(0,2),(1,1),(1,2),(2,2)]
    selected = []; checks = []; moment_checks = 0
    for leg, cut, upper, lower in [(0,0,9,10),(1,3,11,12)]:
        weights = [int(a == leg)+int(b == leg) for a,b in pairs]
        weights += [int(i == leg) for i in range(3)]
        dependent = [s for s in range(9) if roles[s] == 0
                     and any(x*w for x,w in zip(factors[s][:9],weights))]
        assert dependent == ([2,4,5,6,8] if leg == 0 else [2,4,8])
        matching = [(i,s) for i,s in enumerate(sources)
                    if s["id"].startswith(f"lorentz/{leg}/{leg}/")]
        assert len(matching) == 16
        selected.append({"leg": leg, "cut": cut, "upper": upper, "lower": lower,
                         "dependent_ordinary_bases_zero": dependent,
                         "captured_source_rows": [{"ordinal": i, **s} for i,s in matching]})
        for d in [11, 13]:
            for n in range(1,5):
                for s in range(4):
                    base = [0]*16; base[0] = base[3] = 2; base[cut] = n
                    base[upper] = s; base[1] = 3
                    other_upper = 11 if leg == 0 else 9
                    base[other_upper] = 2
                    assert any(all((lo is None or b >= lo) and (hi is None or b <= hi)
                                   for b,(lo,hi) in zip(base,src["domain"])) for _,src in matching)
                    row = defaultdict(Q); row[tuple(base)] += d
                    for slot, factor in enumerate(factors):
                        dc = [x*w for x,w in zip(factor[:9],weights)]
                        coeff = Q(1) if roles[slot] == 2 and base[slot] == 0 else -Q(base[slot])
                        for j in range(9):
                            c = coeff * sum((dc[k]*inv[k][j] for k in range(9)), Q(0))
                            if c:
                                image = base.copy(); image[slot] += 1; image[j] -= 1
                                if any(all((lo is None or v >= lo) and (hi is None or v <= hi)
                                           for v,(lo,hi) in zip(image,z)) for z in zeros): continue
                                row[tuple(image)] += c
                    row = {k:v for k,v in row.items() if v}
                    expected = defaultdict(Q); expected[tuple(base)] = Q(d-2*n)
                    energy = [Q(i == 6+leg) for i in range(9)]
                    for j in range(9):
                        c = Q(-1 if s == 0 else s)*sum((energy[k]*inv[k][j] for k in range(9)),Q(0))
                        if c:
                            image=base.copy();image[upper]+=1;image[j]-=1;expected[tuple(image)]+=c
                    assert row == {k:v for k,v in expected.items() if v}
                    residual = (d-2*n)*moment(d,n,s,0,Q(2,3))
                    residual += (-1 if s==0 else s)*moment(d,n,s+1,1,Q(2,3))
                    assert residual == 0; moment_checks += 1
                    checks.append({"leg":leg,"D":d,"cut":n,"upper":s,
                                   "row": [{"indices":list(k),"coefficient":str(v)} for k,v in sorted(row.items())]})
    result = {"scope":"Exact re-evaluation of the existing production dilation generator using K's saved factor map and original source domains; structural native source descriptors only, not decoded coefficient or rule replay",
              "conclusion":"Already-present local dilation yields D-2n Ward; no additional source policy justified",
              "source_count":len(sources),"generator_row_checks":len(checks),"raised_cut_moment_checks":moment_checks,
              "selected_sources":selected,"checks":checks,
              "inputs_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [META,META.with_suffix('.bin'),ROOT/'src/finite_density/measure.rs',ROOT/'src/finite_density/preparation.rs',Path(__file__)]}}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({"source_count":len(sources),"generator_row_checks":len(checks),"moment_checks":moment_checks,"status":"passed"}))


if __name__ == '__main__': main()
