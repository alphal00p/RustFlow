"""Forcer-type master values for the provided table (pflow-masters.v1) from ACCEPTED route-B samples of the same d0.
A master is defined in the p-integral answers' master_definitions as {"family": F, "index": a, "relation": "M = eps^k * E_G(key)"}
with E_G = I_E / (eps G(1,1))^L. The sample supplies I_E of the integral a (a master value at m = 1 or an --extra target) with its
achieved abs_err; with G11 = eps G(1,1) = Gamma(1+eps) Gamma(1-eps)^2 / Gamma(2-2eps): M = eps^k I_E / G11^L and
err(M) = |eps^k / G11^L| err(I_E) + rounding.
usage: masters_from_sample.py --answers A.json[,B.json] --d0 P/Q --out MASTERS.json --master NAME=SAMPLE.json ...
       [--frame-sha256 NAME=HEX]
Refused (exit 2): a sample of another d0, a rejected or verdict-less sample, a sample whose frame digest differs from the one given,
a target the sample does not contain, a master not in the definitions, an unknown relation."""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pflowio import Refused, MASTERS_SCHEMA, rat, read_answers, read_sample, fmt_err


def main(argv):
    from mpmath import mp, mpf, gamma, nstr, fabs
    mp.dps = 220
    opt = lambda k: argv[argv.index(k) + 1] if k in argv else None
    answers = opt("--answers").split(",")
    d0 = rat(opt("--d0"))
    out = opt("--out")
    frames = dict(x.split("=", 1) for i, x in enumerate(argv) if i > 0 and argv[i - 1] == "--frame-sha256")
    specs = [x for i, x in enumerate(argv) if i > 0 and argv[i - 1] == "--master"]
    if not (answers and out and specs):
        raise Refused(__doc__.split("usage: ")[1].split("\n")[0])
    _, defs, defs_sha, _, sources = read_answers(answers)
    eps_q = (4 - d0) / 2
    eps = mpf(eps_q.numerator) / eps_q.denominator
    G11 = gamma(1 + eps) * gamma(1 - eps) ** 2 / gamma(2 - 2 * eps)
    vals = {}
    for sp in specs:
        name, path = sp.split("=", 1)
        if name not in defs:
            raise Refused(f"master {name} is not defined in the answers")
        md = defs[name]
        m = re.fullmatch(r"M = (?:eps(?:\^(\d+))? \* )?E_G\(key\)", md["relation"].strip())
        if not m:
            raise Refused(f"master {name}: unknown relation {md['relation']!r}")
        k = 0 if md["relation"].strip() == "M = E_G(key)" else int(m.group(1) or 1)
        s = read_sample(path, d0=d0, frame_sha=frames.get(name))
        idx = tuple(md["index"])
        if idx in s["extras"]:
            v, e = s["extras"][idx]
        elif (1, idx) in s["values"]:
            v, e = s["values"][(1, idx)]
        else:
            raise Refused(f"{path}: no value of {md['family']} {list(idx)} (master {name})")
        # E_G divides by (eps G(1,1))^L with L = the loops of the defining family
        L = int(md.get("loops", 0)) or int(re.match(r"P(\d)L", md["family"]).group(1))
        f = eps ** k / G11 ** L  # eps G(1,1) = Gamma(1+eps) Gamma(1-eps)^2 / Gamma(2-2eps) = G11
        M = f * mpf(v)
        err = fabs(f) * mpf(e) + fabs(M) * mpf(10) ** -200
        vals[name] = {"value": nstr(M, 160), "abs_err": fmt_err(err), "sample": s["file"], "sample_sha256": s["sha256"],
                      "frame_sha256": s["frame_sha256"], "family": md["family"], "index": md["index"]}
        sources.append({"file": s["file"], "sha256": s["sha256"], "role": f"sample for master {name}"})
    doc = {"schema": MASTERS_SCHEMA, "d0": f"{d0.numerator}/{d0.denominator}", "master_definitions_sha256": defs_sha,
           "sources": sources, "values": vals, "producer": "RustFlow tools/pflow/masters_from_sample.py"}
    with open(out, "w") as fo:
        json.dump(doc, fo, indent=1)
    for n, x in vals.items():
        print(n, "=", x["value"][:42], "+-", x["abs_err"])


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        sys.exit(2)
