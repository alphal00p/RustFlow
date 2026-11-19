"""Blind PSLQ recognition with independent confirmation (GATE 3).

A dictionary (pflow-pslq-dict.v1) is FROZEN before any value is looked at:
  {"schema": "pflow-pslq-dict.v1", "convention": "MZV indices outer first: zeta(s1,..,sk) = sum_{n1 > .. > nk} n1^-s1 .. nk^-sk",
   "max_height": H, "maxsteps": N,
   "elements": [{"name": "1", "terms": [["1", []]]}, {"name": "z3^2", "terms": [["1", [[3], [3]]]]},
                {"name": "z(3,5)hat", "terms": [["1", [[3, 5]]], ["-29/12", [[8]]]]}, ...]}
  every element is a rational combination of products of MZVs (each MZV an index list).
or generated: hoffman_dictionary(W) = 1 and the Hoffman MZVs zeta(s1..sk), s_i in {2, 3}, of weight <= W (a spanning set of the
MZVs of weight <= W: 20 / 27 / 36 elements through weight 9 / 10 / 11).

Constants: the iterated integral by Hoelder convolution with nested power sums. The two split points 1/2 and 1/3 are ONE algorithm
(they catch truncation errors at the full working precision, not a convention or algebra error). INDEPENDENTLY, every depth-1
factor is compared with mpmath's own zeta and every depth-2 factor with a one-dimensional quadrature
zeta(a,b) = 1/Gamma(a) int_0^1 (-log t)^(a-1) Li_b(t) / (1 - t) dt (mpmath polylog), to INDEPENDENT_DIGITS digits. Depth-3 factors
have no independent evaluation here; the tests check them through stuffle and duality identities.

search(x_A, p1): PSLQ on (x_A, B_1..B_n) ALL rounded to p1 significant digits (20 guard digits), tolerance n H 10^(1 - p1) with H =
max_height (a fixed 10^(3 - p1) rejects every true relation of height above ~10^3); a relation with a nonzero target coefficient
and height <= H is a CANDIDATE. Blind regime: p1 > (n + 1) log10 H (blind_digits): below it relations of height <= H exist by
chance among any n + 1 numbers.
confirm(candidate, x_B, err_B): x_B is an INDEPENDENT value (a disjoint fit, or a higher-precision rerun) with its propagated
absolute uncertainty. ACCEPTED iff |x_B - sum r_b B_b| <= err_B + sum |r_b| delta_b (delta_b: the constants' uncertainty) AND the
independent value carries >= excess_min (default 40) informative digits beyond p1: -log10(err_B/|x_B|) - p1 >= 40.
The search digits are never used for the acceptance. A failed or ambiguous candidate is reported, never dropped.
"""
import hashlib
import json
import sys
from fractions import Fraction as Fr
from functools import lru_cache

from mpmath import mp, mpf, pslq, log10, fabs, nstr

DICT_SCHEMA = "pflow-pslq-dict.v1"


def word_of(ss):
    """MZV index list (outer first) -> iterated-integral word in {0, 1}: 0^(s1-1) 1 0^(s2-1) 1 ..."""
    return "".join("0" * (s - 1) + "1" for s in ss)


def _li_at(word, x, dps):
    """iterated integral I(0; word; x) for 0 < x < 1 by strictly nested power sums (the word ends in 1)"""
    if not word:
        return mpf(1)
    idx, ar = [], 1
    for ch in word:
        if ch == "0":
            ar += 1
        else:
            idx.append(ar)
            ar = 1
    nmax = int(dps / float(-log10(x))) + 40
    partial = [mpf(0) for _ in idx] + [mpf(1)]
    total, scale = mpf(0), mpf(1)
    for n in range(1, nmax):
        scale *= x
        total += scale * partial[1] / mpf(n) ** idx[0]
        for i in range(1, len(idx)):
            partial[i] += partial[i + 1] / mpf(n) ** idx[i]
    return total


@lru_cache(None)
def mzv(ss, dps, split=(1, 2)):
    """zeta(ss) (outer index first) at dps digits, by the Hoelder convolution at p = split:
    I(0; w; 1) = sum_i I(0; dual(prefix_i); 1 - p) I(0; suffix_i; p), dual = reverse + swap 0 <-> 1"""
    with mp.workdps(dps + 20):
        w = word_of(ss)
        p = mpf(split[0]) / split[1]
        tr = str.maketrans("01", "10")
        return +sum(_li_at(w[:i][::-1].translate(tr), 1 - p, dps + 20) * _li_at(w[i:], p, dps + 20) for i in range(len(w) + 1))


INDEPENDENT_DIGITS = 200


@lru_cache(None)
def mzv_independent(ss, dps):
    """zeta(ss) by an algorithm that shares nothing with `mzv`: depth 1 mpmath zeta; depth 2 the polylog quadrature; None for
    depth >= 3"""
    from mpmath import zeta, quad, polylog, log, gamma
    with mp.workdps(dps + 15):
        if len(ss) == 1:
            return +zeta(ss[0])
        if len(ss) == 2:
            a, b = ss
            q = quad(lambda t: (-log(t)) ** (a - 1) * polylog(b, t) / (1 - t), [0, mpf(1) / 4, mpf(1) / 2, mpf(3) / 4, 1])
            return +(q / gamma(a))
    return None


def element_value(el, dps, split=(1, 2)):
    with mp.workdps(dps + 10):
        v = mpf(0)
        for c, prod in el["terms"]:
            t = mpf(Fr(c).numerator) / Fr(c).denominator
            for ss in prod:
                t *= mzv(tuple(ss), dps, split)
            v += t
        return +v


def constants(dic, dps):
    """values of every element (split point 1/2), checked against split point 1/3 at the full precision (truncation) and every
    depth <= 2 factor against `mzv_independent` to min(dps, INDEPENDENT_DIGITS) digits; refuses a disagreement"""
    vals = []
    ind = min(dps, INDEPENDENT_DIGITS)
    for ss in sorted({tuple(f) for el in dic["elements"] for _, prod in el["terms"] for f in prod if len(f) <= 2}):
        a, b = mzv(ss, ind), mzv_independent(ss, ind)
        if fabs(a - b) > fabs(a) * mpf(10) ** -(ind - 5):
            raise ValueError(f"constant zeta{ss}: the independent evaluation differs by {nstr(fabs(a - b), 3)}")
    for el in dic["elements"]:
        a = element_value(el, dps, (1, 2))
        b = element_value(el, dps, (1, 3))
        if a != 0 and fabs(a - b) > fabs(a) * mpf(10) ** -(dps - 5):
            raise ValueError(f"constant {el['name']}: the two evaluations differ by {nstr(fabs(a - b), 3)}")
        vals.append(a)
    return vals


def hoffman_dictionary(W, max_height=10 ** 12, maxsteps=10 ** 6):
    """1 and zeta(s1..sk), s_i in {2,3}, weight <= W"""
    els = [{"name": "1", "terms": [["1", []]]}]

    def rec(ss, w):
        if ss:
            els.append({"name": "z(" + ",".join(map(str, ss)) + ")", "terms": [["1", [list(ss)]]]})
        for s in (2, 3):
            if w + s <= W:
                rec(ss + (s,), w + s)
    rec((), 0)
    return {"schema": DICT_SCHEMA, "convention": "MZV indices outer first", "max_height": max_height, "maxsteps": maxsteps,
            "elements": els, "generated": f"hoffman({W})"}


# algebra generators of the MZVs through weight 11 (one per weight, two at weight 11; prod (1 - x^w)^-g_w = 1/(1 - x^2 - x^3)):
# zeta(2), the odd zetas, and the depth-2/3 irreducibles zeta(5,3), zeta(7,3), zeta(3,5,3) (indices outer first)
STANDARD_GENERATORS = [("z2", [2]), ("z3", [3]), ("z5", [5]), ("z7", [7]), ("z53", [5, 3]), ("z9", [9]), ("z73", [7, 3]), ("z11", [11]),
                       ("z353", [3, 5, 3])]


def standard_dictionary(W, max_height=10 ** 12, maxsteps=10 ** 6):
    """1 and every product of the standard generators of total weight <= W (20 / 27 / 36 elements through 9 / 10 / 11): the
    product basis in which physical coefficients have small heights (unlike the Hoffman words)"""
    gens = [(n, ss, sum(ss)) for n, ss in STANDARD_GENERATORS if sum(ss) <= W]
    els = []

    def rec(start, w, names, prods):
        els.append({"name": "*".join(names) or "1", "terms": [["1", [list(x) for x in prods]]]})
        for i in range(start, len(gens)):
            n, ss, wt = gens[i]
            if w + wt <= W:
                rec(i, w + wt, names + [n], prods + [ss])
    rec(0, 0, [], [])
    return {"schema": DICT_SCHEMA, "convention": "MZV indices outer first", "max_height": max_height, "maxsteps": maxsteps,
            "elements": els, "generated": f"standard({W})"}


def dict_sha(dic):
    return hashlib.sha256(json.dumps(dic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def search_tolerance_log10(p1, max_height, n):
    """a relation of height H among n numbers rounded to p1 digits has its own residual up to ~ n H 10^-p1: the tolerance must
    admit it (a fixed 10^(3 - p1) rejects every true relation of height above ~10^3; measured 18 Nov on GGPPSS int200)"""
    return -p1 + float(log10(mpf(max_height))) + float(log10(mpf(n))) + 1


def blind_digits(n, max_height):
    """search digits above which NO relation of height <= max_height among n + 1 numbers exists by chance (dimension x height;
    Bailey-Borwein): (n + 1) log10(max_height)"""
    return (n + 1) * float(log10(mpf(max_height)))


def search(x, basis_vals, p1, max_height, maxsteps):
    """PSLQ on (x, basis) all rounded to p1 significant digits (the information content), arithmetic with 20 guard digits,
    tolerance search_tolerance_log10, heights <= max_height; returns (relation {index: Fraction} or None, raw integer vector)"""
    n = len(basis_vals) + 1
    with mp.workdps(p1 + 20):
        vec = [mpf(nstr(v, p1)) for v in [x] + list(basis_vals)]
        rel = pslq(vec, tol=mpf(10) ** search_tolerance_log10(p1, max_height, n), maxcoeff=max_height, maxsteps=maxsteps)
    if not rel or rel[0] == 0:
        return None, rel
    return {i: Fr(-r, rel[0]) for i, r in enumerate(rel[1:]) if r != 0}, rel


def confirm(rel, basis_vals, x_b, err_b, p1, const_err, excess_min=40):
    """residual of the candidate on the independent value, its tolerance and the informative excess digits"""
    pred = sum(mpf(c.numerator) / c.denominator * basis_vals[i] for i, c in rel.items())
    res = fabs(x_b - pred)
    tol = mpf(err_b) + sum(fabs(mpf(c.numerator) / c.denominator) for c in rel.values()) * const_err
    inform = float(-log10(mpf(err_b) / fabs(x_b))) if err_b > 0 else 999.0
    excess = inform - p1
    ok = bool(res <= tol and excess >= excess_min)
    return {"residual": nstr(res, 3), "tolerance": nstr(tol, 3), "informative_digits": inform, "excess_digits": excess,
            "accepted": ok, "reason": "ok" if ok else ("residual above the propagated uncertainty" if res > tol else f"excess {excess:.1f} < {excess_min}")}


def recognize(x_a, x_b, err_b, dic, p1, excess_min=40, const_dps=None):
    """the full protocol: search on x_A at p1, confirm on x_B; also re-searches at p1 + 10 and p1 + 20 (ambiguity report)"""
    dps = const_dps or int(-float(log10(mpf(err_b) / fabs(x_b)))) + 30
    mp.dps = dps
    vals = constants(dic, dps)
    const_err = mpf(10) ** -(dps - 10)
    rel, raw = search(x_a, vals, p1, dic["max_height"], dic["maxsteps"])
    out = {"dictionary_sha256": dict_sha(dic), "elements": [e["name"] for e in dic["elements"]], "p1": p1}
    if rel is None:
        out.update(verdict="NO CANDIDATE", raw=[int(r) for r in raw] if raw else None)
        return out
    out["candidate"] = {dic["elements"][i]["name"]: str(c) for i, c in rel.items()}
    out["height"] = max(abs(int(r)) for r in raw)
    out["confirmation"] = confirm(rel, vals, x_b, err_b, p1, const_err, excess_min)
    alt = []
    for q in (p1 + 10, p1 + 20):
        r2, _ = search(x_a, vals, q, dic["max_height"], dic["maxsteps"]) if q <= dps - 10 else (None, None)
        if r2 is not None and r2 != rel:
            alt.append({"p": q, "candidate": {dic["elements"][i]["name"]: str(c) for i, c in r2.items()}})
    out["ambiguous_at"] = alt
    out["verdict"] = "REJECTED" if not out["confirmation"]["accepted"] else ("AMBIGUOUS" if alt else "ACCEPTED")
    return out


def scan(x, dic, digits, truth=None, const_dps=None):
    """discovery scan: at each p the candidate and whether it equals `truth` ({name: Fraction}) or, without truth, whether it is
    confirmed by x itself at full precision (x must then carry >= max(digits) + 40 digits)"""
    dps = const_dps or max(digits) + 60
    mp.dps = dps
    vals = constants(dic, dps)
    rows = []
    for p in digits:
        rel, raw = search(x, vals, p, dic["max_height"], dic["maxsteps"])
        if rel is None:
            rows.append({"p": p, "status": "no target relation"})
            continue
        named = {dic["elements"][i]["name"]: c for i, c in rel.items()}
        if truth is not None:
            st = "exact" if named == truth else "other"
        else:
            pred = sum(mpf(c.numerator) / c.denominator * vals[i] for i, c in rel.items())
            st = "confirmed" if fabs(x - pred) < fabs(x) * mpf(10) ** -(p + 40) else "other"
        rows.append({"p": p, "status": st, "height": max(abs(int(r)) for r in raw), "support": len(rel)})
    return rows


USAGE = """usage: recognize.py (--dict DICT.json --dict-sha256 HEX | --standard W | --hoffman W) [--dict-sha256 HEX] [--max-height H]
         --p1 P [--excess 40]
         (--fit-a FIT_A.json --fit-b FIT_B.json --coef K | --x-a VALUE --x-b VALUE --err-b ABS_ERR) [--out RESULT.json]
FIT_A / FIT_B: SEALED pflow-fit.v1 outputs (fitgrid.py --out) of the same coefficient (same target, frame, scheme, loops, pole
order) on DISJOINT node sets (no d0 and no sample in common); fit A must carry >= P informative digits. The search uses fit A's
coefficient c_K rounded to P digits, the confirmation fit B's c_K with its propagated abs_err. Prints the verdict (ACCEPTED / REJECTED / AMBIGUOUS /
NO CANDIDATE) and writes the full record (dictionary sha256, candidate, residual, tolerance, excess digits)."""


def main(argv):
    opt = lambda k: argv[argv.index(k) + 1] if k in argv else None
    H = int(float(opt("--max-height") or "1e12"))
    if opt("--dict"):
        dic = json.load(open(opt("--dict")))
        if dic.get("schema") != DICT_SCHEMA:
            raise SystemExit(f"REFUSED: {opt('--dict')}: schema {dic.get('schema')!r}")
    elif opt("--standard"):
        dic = standard_dictionary(int(opt("--standard")), H)
    elif opt("--hoffman"):
        dic = hoffman_dictionary(int(opt("--hoffman")), H)
    else:
        raise SystemExit(USAGE)
    # the dictionary must be pinned: committed with its sha256 before any value is looked at (required for a dictionary file)
    pin = opt("--dict-sha256")
    if opt("--dict") and not pin:
        raise SystemExit("REFUSED: --dict needs --dict-sha256 (the pinned digest of the frozen dictionary)")
    if pin and pin != dict_sha(dic):
        raise SystemExit(f"REFUSED: dictionary sha256 {dict_sha(dic)[:16]}... is not the pinned {pin[:16]}...")
    p1 = int(opt("--p1") or 0)
    if p1 <= 0:
        raise SystemExit(USAGE)
    mp.dps = 2 * p1 + 100
    if opt("--fit-a"):
        sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
        from pflowio import check_seal, Refused
        fa, fb = json.load(open(opt("--fit-a"))), json.load(open(opt("--fit-b")))
        try:
            check_seal(opt("--fit-a"), fa)
            check_seal(opt("--fit-b"), fb)
        except Refused as e:
            raise SystemExit(f"REFUSED: {e}")
        k = opt("--coef")
        sa, sb = fa.get("spec", {}), fb.get("spec", {})
        for f in ("target", "frame_sha256", "scheme", "loops", "pole_order"):
            if sa.get(f) != sb.get(f):
                raise SystemExit(f"REFUSED: fits A and B differ in {f} ({sa.get(f)!r} vs {sb.get(f)!r}): not the same coefficient")
        if set(sa.get("nodes", {})) & set(sb.get("nodes", {})):
            raise SystemExit("REFUSED: fits A and B share nodes (d0): the confirmation would not be independent")
        used_a = {u["sha256"] for u in fa.get("used_samples", [])}
        used_b = {u["sha256"] for u in fb.get("used_samples", [])}
        if used_a & used_b:
            raise SystemExit("REFUSED: fit A and fit B share samples: the confirmation would not be independent")
        x_a, x_b, err_b = mpf(fa["coefficients"][k]["value"]), mpf(fb["coefficients"][k]["value"]), mpf(fb["coefficients"][k]["abs_err"])
        err_a = mpf(fa["coefficients"][k]["abs_err"])
        if err_a > 0 and x_a != 0 and float(-log10(err_a / fabs(x_a))) < p1:
            raise SystemExit(f"REFUSED: fit A carries {float(-log10(err_a / fabs(x_a))):.1f} informative digits < the search digits {p1}")
    elif opt("--x-a"):
        x_a, x_b, err_b = mpf(opt("--x-a")), mpf(opt("--x-b")), mpf(opt("--err-b"))
    else:
        raise SystemExit(USAGE)
    res = recognize(x_a, x_b, err_b, dic, p1, excess_min=float(opt("--excess") or 40))
    res["blind_digits"] = blind_digits(len(dic["elements"]), dic["max_height"])
    print(res["verdict"], json.dumps({k: res[k] for k in ("candidate", "confirmation", "ambiguous_at") if k in res}))
    if opt("--out"):
        json.dump(res, open(opt("--out"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main(sys.argv[1:])
