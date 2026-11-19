"""Laurent coefficients of one target from a grid of AUTHENTICATED route-B samples (pflow-fit.v1).
usage: fitgrid.py SPEC.json SAMPLE.json ... [--out FIT.json]
SPEC (pflow-fit-spec.v1):
  {"schema": "pflow-fit-spec.v1", "nodes": {"p/q": SHA256, ...}   (every node with the sha256 of its ACCEPTED sample file; ALL
   required, each exactly once), "frame_sha256": HEX (required), "need_target_digits": T, "need_abs_err": X,
   "target": {"name": N, "master": [a...], "m": 1} | {"name": N, "extra": [a...]}, "loops": L, "scheme": "G" | "E",
   "pole_order": p, "production": true | false (REQUIRED; nothing else is accepted), "legacy_samples": {"p/q": SHA256} (only the 12
   archived 4L nodes; loops = 4)}
The contract class (production or not) is printed on every result line.
Every sample is authenticated (pflowio.read_sample): content seal, the file's sha256 equal to the pinned one, exact d0 in the node
list, the frame digest recomputed and equal to frame_sha256, ACCEPTED with a verdict, effective target >= need_target_digits, the
target's achieved abs_err (after the normalisation to the scheme, times |eps|^p) <= need_abs_err. A production contract refuses
every downgraded sample (no two-precision pass, unreachable inputs, legacy). Missing, duplicate, extra or unpinned nodes refuse
the fit (exit 2, "REFUSED: ..."), never a smaller fit. legacy_samples is refused at loops >= 5 and for any file that is not one of
the 12 archived four-loop nodes; every printed line of a fit that used one says LEGACY.
The fit: y_i = eps_i^p V(eps_i) is interpolated by the degree n-1 polynomial through all n nodes, c_k = sum_i w_ik y_i with the
EXACT rational cardinal weights; abs_err(c_k) = sum_i |w_ik| err(y_i) (worst case, correlations not assumed) and
amplification_log10 = log10 sum_i |w_ik| (from the node set and p alone). The output records the digests of the samples used."""
import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pflowio import Refused, rat, read_sample, digest, sha256_file, fmt_err, load_json, LEGACY_4L_NODES, seal

FIT_SPEC_SCHEMA = "pflow-fit-spec.v1"


def cardinal_weights(nodes, k):
    """w_i = [x^k] l_i(x), l_i the Lagrange cardinal polynomials of the exact rational nodes"""
    ws = []
    for i, x in enumerate(nodes):
        poly = [Fr(1)]
        for j, b in enumerate(nodes):
            if i != j:
                new = [Fr(0)] * (len(poly) + 1)
                for h, v in enumerate(poly):
                    new[h] += -b * v / (x - b)
                    new[h + 1] += v / (x - b)
                poly = new
        ws.append(poly[k] if k < len(poly) else Fr(0))
    return ws


def g11(eps):
    from mpmath import gamma
    return gamma(1 + eps) * gamma(1 - eps) ** 2 / gamma(2 - 2 * eps)


def fit(spec, files, mp_dps=None):
    from mpmath import mp, mpf, fabs, log10, nstr
    if spec.get("schema") != FIT_SPEC_SCHEMA:
        raise Refused(f"fit spec schema {spec.get('schema')!r}, expected {FIT_SPEC_SCHEMA!r}")
    if "allow_legacy" in spec:
        raise Refused("fit spec: allow_legacy was replaced by legacy_samples {d0: sha256} (the 12 archived four-loop nodes only)")
    if not isinstance(spec.get("nodes"), dict):
        raise Refused("fit spec: nodes must pin every node's sample: {\"p/q\": sha256, ...}")
    pins = {rat(k): v for k, v in spec["nodes"].items()}
    if len(pins) != len(spec["nodes"]):
        raise Refused("fit spec: duplicate node")
    nodes = list(pins)
    if not isinstance(spec.get("frame_sha256"), str):
        raise Refused("fit spec: frame_sha256 is required")
    production = spec.get("production")
    if not isinstance(production, bool):
        raise Refused(f"fit spec: production must be true or false (got {production!r}): the contract class is never implicit")
    legacy = {rat(k): v for k, v in (spec.get("legacy_samples") or {}).items()}
    if legacy:
        if int(spec["loops"]) != 4:
            raise Refused(f"fit spec: legacy_samples at loops = {spec['loops']}: the legacy path exists only for the 12 archived four-loop nodes")
        for k, v in legacy.items():
            if v not in LEGACY_4L_NODES or pins.get(k) != v:
                raise Refused(f"fit spec: legacy sample {k} ({v[:16]}...) is not one of the 12 archived four-loop nodes pinned in nodes")
        if production:
            raise Refused("fit spec: a production contract cannot use legacy_samples")
    tgt = spec["target"]
    key = ("extra", tuple(tgt["extra"])) if "extra" in tgt else ("value", (tgt.get("m", 1), tuple(tgt["master"])))
    p, L, scheme = int(spec["pole_order"]), int(spec["loops"]), spec.get("scheme", "G")
    need_t, need_e = spec.get("need_target_digits"), spec.get("need_abs_err")
    if need_t is None or need_e is None:
        raise Refused("fit spec: need_target_digits and need_abs_err are required")
    got = {}
    for f in files:
        s = read_sample(f, frame_sha=spec["frame_sha256"], need_target=float(need_t), legacy_pins=set(legacy.values()))
        if s["d0"] not in nodes:
            raise Refused(f"{f}: d0 = {s['d0']} is not a node of the fit spec")
        if s["sha256"] != pins[s["d0"]]:
            raise Refused(f"{f}: sha256 {s['sha256'][:16]}... is not the sample pinned for d0 = {s['d0']} ({pins[s['d0']][:16]}...)")
        if production and s["downgrades"]:
            raise Refused(f"{f}: downgraded sample refused by a production contract: {'; '.join(s['downgrades'])}")
        if s["d0"] in got:
            raise Refused(f"{f}: a second sample for d0 = {s['d0']} ({got[s['d0']]['file']})")
        got[s["d0"]] = s
    missing = [str(n) for n in nodes if n not in got]
    if missing:
        raise Refused(f"missing nodes {missing}: the spec requires all {len(nodes)} (no smaller fit)")
    mp.dps = mp_dps or 60 + int(max(float(need_t), 40)) * 2
    eps_q = [(4 - n) / 2 for n in nodes]
    ys, es, used = [], [], []
    for n, e in zip(nodes, eps_q):
        s = got[n]
        tab = s["extras"] if key[0] == "extra" else s["values"]
        if key[1] not in tab:
            raise Refused(f"{s['file']}: no value for the target {tgt['name']}")
        v, ae = tab[key[1]]
        ef = mpf(e.numerator) / e.denominator
        norm = 1 / g11(ef) ** L if scheme == "G" else mpf(1)
        y = ef ** p * mpf(v) * norm
        ey = fabs(ef ** p * norm) * mpf(ae)
        if ey > mpf(need_e):
            raise Refused(f"{s['file']}: target uncertainty {nstr(ey, 3)} (normalised) above the need {need_e}")
        ys.append(y)
        es.append(ey)
        used.append({"file": s["file"], "sha256": s["sha256"], "d0": f"{n.numerator}/{n.denominator}", "legacy": s["legacy"],
                     "target_digits": s["target_digits"], "downgrades": s["downgrades"]})
    out = {}
    for k in range(-p, len(nodes) - p):
        w = cardinal_weights(eps_q, k + p)
        c = sum(mpf(wi.numerator) / wi.denominator * y for wi, y in zip(w, ys))
        amp = sum(abs(wi) for wi in w)
        err = sum(fabs(mpf(wi.numerator) / wi.denominator) * ey for wi, ey in zip(w, es))
        out[k] = {"value": nstr(c, mp.dps - 10), "abs_err": fmt_err(err), "amplification_log10": float(log10(mpf(amp.numerator) / amp.denominator)),
                  "c": c, "err": err}
    return out, used


def main(argv):
    from mpmath import nstr
    opt = lambda k: argv[argv.index(k) + 1] if k in argv else None
    pos = [a for i, a in enumerate(argv) if not a.startswith("--") and (i == 0 or argv[i - 1] != "--out")]
    spec = load_json(pos[0])
    cs, used = fit(spec, pos[1:])
    flag = "  [PRODUCTION contract]" if spec["production"] else "  [NOT a production contract]"
    flag += "  LEGACY (frame unverified)" if any(u["legacy"] for u in used) else ""
    flag += "  DOWNGRADED inputs" if any(u["downgrades"] for u in used) and "LEGACY" not in flag else ""
    for k, c in cs.items():
        print(f"c_{k:+d} = {nstr(c['c'], 30)}  +- {c['abs_err']}  (amplification 10^{c['amplification_log10']:.3f}){flag}")
    if opt("--out"):
        doc = {"schema": "pflow-fit.v1", "spec": spec, "spec_sha256": digest(spec), "used_samples": used,
               "coefficients": {str(k): {x: c[x] for x in ("value", "abs_err", "amplification_log10")} for k, c in cs.items()}}
        json.dump(seal(doc), open(opt("--out"), "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        sys.exit(2)
