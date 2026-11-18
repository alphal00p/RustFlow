"""p-integral's d-exact reductions (num-support.req-answer.v1) -> the provided table (pflow-provided.v2) of the RustFlow propagator
mode at ONE exact rational d0 (E-normalisation: Euclidean, measure d^dk/pi^(d/2), q^2 = 1).
    (-1)^(sum a) E_G(I) = sum_terms coefficient(d) * basis,  E_G = I_E / (eps G(1,1))^L,  basis = G[(m,e),...] = prod Gamma(1+m eps)^e
    or M:<name> (a master as defined in the answers' master_definitions; values from MASTERS.json).
usage: answer2prov.py KEYS.json D0 OUT.json ANSWER.json [ANSWER.json ...] [--masters MASTERS.json] [--dps 160] [--need-digits N]

Authentication (every failure exits 2 with "REFUSED: ..." and writes nothing):
- every answer file has the same master definitions and meaning (their digest goes into the table);
- the masters file names the SAME exact d0 (and, for pflow-masters.v1, the same definitions digest);
- coefficients are evaluated EXACTLY in rationals at the exact d0 and cross-checked against the answer's own exact coefficients at
  that d0 when it lists them;
- every row carries its ACHIEVED absolute uncertainty: |G11^L| (sum |c_t| |G_t| err(M_t) + evaluation rounding) + printing
  rounding; --need-digits N refuses a table whose worst row has fewer than N relative digits.
Rows built from a legacy master file (no measured uncertainty) are listed in "unverified"; RustFlow refuses such a table unless
the run explicitly allows unverified provided values."""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pflowio import (Refused, PROVIDED_SCHEMA, NORMALISATION, rat, eval_exact, read_answers, read_masters, fmt_err, rel_digits,
                     sha256_file)


def main(argv):
    from mpmath import mp, mpf, gamma, nstr, fabs
    if "--master-digits" in argv:
        raise Refused("--master-digits was removed: master uncertainties come from the masters file and are propagated")
    opt = lambda k, d=None: argv[argv.index(k) + 1] if k in argv else d
    flags = ("--dps", "--masters", "--need-digits")
    pos = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or argv[i - 1] not in flags)]
    if len(pos) < 4:
        raise Refused(__doc__.split("\n")[3])
    keys, d0s, out, answers = pos[0], pos[1], pos[2], pos[3:]
    dps = int(opt("--dps", "160"))
    need = opt("--need-digits")
    mfile = opt("--masters")
    mp.dps = dps + 20
    d0 = rat(d0s)
    eps_q = (4 - d0) / 2
    eps = mpf(eps_q.numerator) / eps_q.denominator
    ans, defs, defs_sha, meaning_sha, sources = read_answers(answers)
    masters, verified = {}, {}
    if mfile:
        mv, msrc = read_masters(mfile, d0, defs_sha)
        sources += msrc
        for k, (v, e, ok) in mv.items():
            if k not in defs:
                raise Refused(f"{mfile}: master {k} is not defined in the answers' master_definitions")
            masters[k] = (mpf(v), mpf(e))
            verified[k] = ok
    G11 = gamma(1 + eps) * gamma(1 - eps) ** 2 / gamma(2 - 2 * eps)

    def gtag(s):
        v = mpf(1)
        for m, e in re.findall(r"\((-?\d+),(-?\d+)\)", s):
            v *= gamma(1 + int(m) * eps) ** int(e)
        return v

    K = json.load(open(keys))
    id2key = {v: k for k, v in K["ids"].items()}
    table, unverified, n_m, worst = [], [], 0, None
    for f, A in ans:
        for x in A["integrals"]:
            if x["id"] not in id2key:
                continue
            # the answer's own exact coefficients at this d0, if listed: an independent check of the evaluation
            at = next((v for k, v in (x.get("at_d0") or {}).items() if rat(k) == d0), None)
            tot, err, mag, uses_unverified = mpf(0), mpf(0), mpf(0), False
            for i, t in enumerate(x["terms"]):
                b = t["basis"]
                c = eval_exact(t["coefficient"], {"d": d0})
                if at is not None and rat(at["terms"][i]["coefficient_exact"]) != c:
                    raise Refused(f"{f}: {x['id']} term {i}: coefficient at d0 = {c}, the answer lists {at['terms'][i]['coefficient_exact']}")
                cf = mpf(c.numerator) / c.denominator
                gv = gtag(b[b.index("G["):]) if "G[" in b else mpf(1)
                mm = re.match(r"M:(\w+)", b)
                if mm:
                    name = mm.group(1)
                    if name not in masters:
                        raise Refused(f"value of master {name} needed (--masters)")
                    mv, me = masters[name]
                    tot += cf * mv * gv
                    err += fabs(cf) * fabs(gv) * me
                    mag += fabs(cf * mv * gv)
                    uses_unverified |= not verified[name]
                    n_m += 1
                else:
                    tot += cf * gv
                    mag += fabs(cf * gv)
            key = json.loads(id2key[x["id"]])
            L = len(key) - len({v for e in key for v in e[:2]}) + 1
            s = fabs(G11) ** L
            val = x["sign"] * tot * G11 ** L
            ae = s * (err + mag * mpf(10) ** -(dps + 15)) + fabs(val) * mpf(10) ** -(dps - 1)
            dg = rel_digits(val, ae)
            worst = dg if worst is None else min(worst, dg)
            if uses_unverified:
                unverified.append(x["id"])
            table.append({"edges": [[e[0], e[1], e[2]] for e in key], "ext": [0, 1], "value": nstr(val, dps), "abs_err": fmt_err(ae),
                          "id": x["id"], "family": x.get("family"), "index": x.get("index"), "source": Path(f).name})
    if need is not None and worst is not None and worst < float(need):
        raise Refused(f"the worst row has {worst:.1f} relative digits < --need-digits {need}")
    doc = {"schema": PROVIDED_SCHEMA, "normalisation": NORMALISATION, "d0": f"{d0.numerator}/{d0.denominator}",
           "master_definitions_sha256": defs_sha, "answer_meaning_sha256": meaning_sha, "sources": sources,
           "keys": {"file": str(keys), "sha256": sha256_file(keys)}, "unverified": unverified,
           "worst_relative_digits": worst, "producer": "RustFlow tools/pflow/answer2prov.py", "values": table}
    with open(out, "w") as fo:
        json.dump(doc, fo, indent=1)
    print(len(table), "provided values at d0", d0, f"({n_m} master terms; worst {worst:.1f} relative digits; {len(unverified)} unverified rows)")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        sys.exit(2)
