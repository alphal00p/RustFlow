"""Provenance and input authentication for the pflow (route B) pipeline: digests, exact rational arithmetic at the exact d0,
and the checked readers of every file that crosses a stage boundary (p-integral answers, master values, provided tables,
route-B samples). Every reader REFUSES (raises Refused) instead of guessing: a d0 mismatch, a digest mismatch, a rejected or
verdict-less sample, a missing uncertainty or one above the need.

Uncertainties are ABSOLUTE and achieved (propagated), never a nominal digit stamp:
- a provided value sum_i c_i B_i carries sum_i |c_i| |B_i| err_i plus the evaluation and printing rounding;
- a route-B sample value carries 10^-digits_estimate max(|v|, 10^-30 scale) (digits_estimate = min over every gate quantity)
  or the two-precision difference, whichever is larger (pflow-route-b.v2 writes it as `abs_err`).
"""
import ast
import hashlib
import json
import math
from fractions import Fraction as Fr
from pathlib import Path

PROVIDED_SCHEMA = "pflow-provided.v2"
MASTERS_SCHEMA = "pflow-masters.v1"
SAMPLE_SCHEMA = "pflow-route-b.v2"
ANSWER_SCHEMA = "num-support.req-answer.v1"
NORMALISATION = "E, P^2=1"


# the ONLY route-B v1 samples a fit may ever accept (allow_legacy): the 12 archived four-loop nodes of the blind eps^1 test
# (s1/runs/rf4/grid-r, sha256 of the files). Any other v1 file - in particular any five-loop sample - is refused.
LEGACY_4L_NODES = {
    "c9a17c5627ddb3f2a4ce7bef5b8ac4f910cec875acf51c035446af588ce3f90c": "node-1003967_250000.json",
    "205d5f33181240745a3add3640dbac6cc6f049edd581e7d4923febb149e669e2": "node-1004957_250000.json",
    "b07cbb798ca34f0b65d04cf2e9bcd47389dcb1fefaf3a88b7a4c624f36e64e0e": "node-1990761_500000.json",
    "fb8c3251b688be552e397b1b4bb1e931138c4c18758cf5e3e88265744bc46611": "node-1996173_500000.json",
    "a5cebf6ee02660ef0c502ada8e7a6edc1bc7fcc780005839b16779bcea1af968": "node-2003827_500000.json",
    "0fc3f6de6f87f4eb0d9ccde3020eb64570cff1fa8a93820688784050512327a9": "node-2009239_500000.json",
    "1a377ba8413f91c3c198dddfa6c8e97af2e3fbb4a5be0642e5c080d8cfbd3926": "node-249239_62500.json",
    "38c3d767186396c378b00b4e3e585eb1fb551097e11fc09beb259325f2ae42ac": "node-250761_62500.json",
    "0e4895401b0c17243ff5ad7be81720431df2356cb23306a092ad27a2af558c8d": "node-399739_100000.json",
    "7bcd5b6b5cf4b77b14f0e9054c8293a02843cc6335ac6c3c4344753ee5cd9377": "node-400261_100000.json",
    "485d2a01aa1261485bc0eb3a28203210de822fc982202102e7736657ce2ede4f": "node-995043_250000.json",
    "0225aff5c46a4690e3cd7c83724d81f69812d992136faa92392262da932ac162": "node-996033_250000.json",
}


class Refused(Exception):
    """an input failed authentication; the message names the file and the check"""


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def digest(obj):
    """sha256 of the canonical JSON of obj (sorted keys, no spaces)"""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def rat(s):
    """exact rational from "p/q", an int, or a Fraction; floats are refused"""
    if isinstance(s, Fr):
        return s
    if isinstance(s, bool) or isinstance(s, float):
        raise Refused(f"not an exact rational: {s!r}")
    if isinstance(s, int):
        return Fr(s)
    if isinstance(s, str) and all(c in "0123456789-/ " for c in s) and s.strip():
        return Fr(s.replace(" ", ""))
    raise Refused(f"not an exact rational: {s!r}")


def eval_exact(expr, env):
    """exact value of an arithmetic expression in integers and the names of env (Fractions): + - * / ^ (integer powers), parentheses.
    Every node is evaluated in Fractions (no binary floating point anywhere)."""
    tree = ast.parse(expr.replace("^", "**"), mode="eval")

    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, int) and not isinstance(n.value, bool):
            return Fr(n.value)
        if isinstance(n, ast.Name) and n.id in env:
            return env[n.id]
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = ev(n.operand)
            return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                if b == 0:
                    raise Refused(f"division by zero in {expr!r} at {env}")
                return a / b
            if isinstance(n.op, ast.Pow):
                if b.denominator != 1:
                    raise Refused(f"non-integer power in {expr!r}")
                return a ** int(b)
        raise Refused(f"unsupported element {ast.dump(n)[:60]} in {expr!r}")

    return ev(tree)


def rel_digits(value, abs_err):
    """-log10(abs_err / |value|) (inf for an exact value, 0 for a zero value with an error)"""
    from mpmath import mpf, log10, fabs
    v, e = fabs(mpf(value)), mpf(abs_err)
    if e == 0:
        return math.inf
    if v == 0:
        return 0.0
    return float(-log10(e / v))


def fmt_err(x):
    """an uncertainty printed with 3 significant digits, rounded UP"""
    from mpmath import mpf, nstr, log10, floor, ceil
    x = mpf(x)
    if x == 0:
        return "0"
    k = int(floor(log10(x))) - 2
    m = ceil(x / mpf(10) ** k)
    return nstr(m * mpf(10) ** k, 3)


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError) as e:
        raise Refused(f"{path}: {e}")


def check_d0(path, got, want):
    if rat(got) != rat(want):
        raise Refused(f"{path}: d0 = {got}, required d0 = {want} (exact)")


# ---------------- p-integral answers ----------------
def read_answers(paths):
    """the answer files of one request: same master definitions and meaning in every file. Returns (answers, definitions,
    definitions_sha256, meaning_sha256, sources)"""
    out, defs, meaning, sources = [], None, None, []
    for p in paths:
        a = load_json(p)
        if a.get("schema") != ANSWER_SCHEMA:
            raise Refused(f"{p}: schema {a.get('schema')!r}, expected {ANSWER_SCHEMA!r}")
        md = a.get("master_definitions")
        if not isinstance(md, dict):
            raise Refused(f"{p}: master_definitions missing")
        if defs is None:
            defs, meaning = md, a.get("meaning")
        elif digest(md) != digest(defs) or a.get("meaning") != meaning:
            raise Refused(f"{p}: master definitions / meaning differ from {paths[0]}")
        out.append((p, a))
        sources.append({"file": str(p), "sha256": sha256_file(p), "role": "answer"})
    return out, defs, digest(defs), digest(meaning), sources


# ---------------- master values ----------------
def read_masters(path, d0, defs_sha=None):
    """{name: (value str, abs_err str, verified bool)} from a pflow-masters.v1 file (d0, definitions digest and sources checked)
    or from a LEGACY file {name: value, "d0", "from"}: its d0 must match; its values have no measured uncertainty and are marked
    unverified (abs_err = |v| 10^-(w - 20) of the source run, if the source exists; refused otherwise)."""
    m = load_json(path)
    if "d0" not in m:
        raise Refused(f"{path}: no d0: a master file must name its dimension")
    check_d0(path, m["d0"], d0)
    src = [{"file": str(path), "sha256": sha256_file(path), "role": "masters"}]
    if m.get("schema") == MASTERS_SCHEMA:
        if defs_sha is not None and m.get("master_definitions_sha256") != defs_sha:
            raise Refused(f"{path}: master_definitions_sha256 {m.get('master_definitions_sha256')} != the answers' {defs_sha}")
        vals = {}
        for k, x in m["values"].items():
            if "abs_err" not in x:
                raise Refused(f"{path}: master {k} has no abs_err")
            vals[k] = (x["value"], x["abs_err"], True)
        return vals, src + m.get("sources", [])
    if m.get("schema") is not None:
        raise Refused(f"{path}: unknown schema {m.get('schema')!r}")
    # legacy file
    from mpmath import mpf, fabs
    w = None
    fr = m.get("from")
    if fr and Path(fr).exists():
        s = load_json(fr)
        check_d0(fr, s.get("d0"), d0)
        w = s.get("config", {}).get("w")
        src.append({"file": str(fr), "sha256": sha256_file(fr), "role": "masters source (legacy, no verdict)"})
    if w is None:
        raise Refused(f"{path}: legacy master file without a readable source run ('from'): no uncertainty can be assigned")
    vals = {}
    for k, v in m.items():
        if k in ("d0", "from"):
            continue
        vals[k] = (v, fmt_err(fabs(mpf(v)) * mpf(10) ** -(int(w) - 20)), False)
    return vals, src


# ---------------- route-B samples ----------------
def read_sample(path, d0=None, frame_sha=None, need_target=None, allow_legacy=False):
    """a route-B sample, authenticated: schema, exact d0 (if given), frame digest (if given), ACCEPTED with a verdict, the
    effective target >= need_target. Returns dict(d0, frame_sha256, target_digits, values: {(m, tuple(master)): (value, abs_err)},
    extras: {tuple(target): (value, abs_err)}, sha256, legacy).
    allow_legacy: a pflow-route-b.v1 sample (no frame, no abs_err) is accepted only if it is one of LEGACY_4L_NODES (by sha256) and
    has accepted == true; its uncertainty is
    the bound its passed gates certify, |v| 10^-T (T = min(target or w - gate_margin, digits_estimate if present)), and its frame
    is unverified (legacy = True)."""
    s = load_json(path)
    sch = s.get("schema")
    if d0 is not None:
        check_d0(path, s.get("d0"), d0)
    if "accepted" not in s:
        raise Refused(f"{path}: no acceptance verdict")
    if s["accepted"] is not True:
        fails = (s.get("acceptance") or s.get("diag", {}).get("acceptance", {})).get("failures")
        raise Refused(f"{path}: REJECTED sample (accepted = {s['accepted']}; failures {fails})")
    out = {"file": str(path), "sha256": sha256_file(path), "d0": rat(s["d0"]), "values": {}, "extras": {}}
    if sch == SAMPLE_SCHEMA:
        if s.get("accept_diagnostics"):
            raise Refused(f"{path}: accepted only through --accept-diagnostics")
        out["frame_sha256"] = s["frame_sha256"]
        out["target_digits"] = float(s["target_digits"])
        out["legacy"] = False
        for v in s["values"]:
            if "abs_err" not in v:
                raise Refused(f"{path}: value without abs_err")
            out["values"][(v["m"], tuple(v["master"]))] = (v["value"], v["abs_err"])
        for x in s.get("extras", []):
            if "error" in x:
                continue
            out["extras"][tuple(x["target"])] = (x["value"], x["abs_err"])
    elif sch == "pflow-route-b.v1" and allow_legacy:
        from mpmath import mpf, fabs
        if out["sha256"] not in LEGACY_4L_NODES:
            raise Refused(f"{path}: a legacy v1 sample is accepted only for the 12 archived four-loop nodes (s1/runs/rf4/grid-r); this file is not one of them")
        c = s.get("config", {})
        t = c.get("target_digits")
        T = float(t) if t is not None else float(c["w"]) - float(c.get("gate_margin", 20))
        if s.get("digits_estimate") is not None:
            T = min(T, float(s["digits_estimate"]))
        out["frame_sha256"] = None
        out["target_digits"] = T
        out["legacy"] = True
        for v in s["values"]:
            out["values"][(v["m"], tuple(v["master"]))] = (v["value"], fmt_err(fabs(mpf(v["value"])) * mpf(10) ** -T))
    else:
        raise Refused(f"{path}: schema {sch!r} (expected {SAMPLE_SCHEMA!r}{' or a legacy v1' if allow_legacy else ''})")
    if frame_sha is not None and out["frame_sha256"] != frame_sha:
        raise Refused(f"{path}: frame_sha256 {out['frame_sha256']} != required {frame_sha}")
    if need_target is not None and out["target_digits"] < need_target:
        raise Refused(f"{path}: accepted at target {out['target_digits']:.1f} digits < the {need_target} needed")
    return out
