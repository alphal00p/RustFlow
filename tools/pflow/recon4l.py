"""four-loop eps^1 blind test (no1 / no6, G-scheme) from the 12-node grid, through the authenticated fit (fitgrid.py) and the
contracts specs/p4l-<name>-blind-eps1.json: every node must be present once, pinned by sha256, accepted, at target >= 18 digits; a
missing, duplicate, extra, unpinned or REJECTED sample refuses the fit (exit 2). These 12 archived nodes are LEGACY v1 samples (the
only ones the legacy path admits); every printed line says so. c_-1 and c_0 are compared with the exact expansions; c_1 is BLIND.
The sub-grid value (the two outermost nodes dropped) is printed as a self-consistency diagnostic only.
usage: recon4l.py GRIDDIR [--out PREFIX]"""
import glob
import json
import sys
from pathlib import Path

H = Path(__file__).resolve().parent
sys.path.insert(0, str(H))
from pflowio import Refused, load_json, seal
from fitgrid import fit


def main(argv):
    from mpmath import mp, mpf, zeta, nstr, log10, fabs
    files = sorted(glob.glob(argv[0] + "/node-*.json"))
    z = zeta
    exact = {"no1": {-1: lambda: -5 * z(5), 0: lambda: 45 * z(5) - 41 * z(3) ** 2 - mpf(25) / 2 * z(6) + mpf(161) / 2 * z(7)},
             "no6": {-1: lambda: -5 * z(5), 0: lambda: 45 * z(5) - 17 * z(3) ** 2 - mpf(25) / 2 * z(6)}}
    for name in ("no1", "no6"):
        spec = load_json(H / "specs" / f"p4l-{name}-blind-eps1.json")
        cs, used = fit(spec, files, mp_dps=120)
        # sub-grid diagnostic: the two outermost |eps| nodes dropped
        sub = dict(spec)
        order = sorted(spec["nodes"], key=lambda d: abs(4 - eval(d)))
        sub["nodes"] = {k: spec["nodes"][k] for k in order[:-2]}
        if spec.get("legacy_samples"):
            sub["legacy_samples"] = {k: v for k, v in spec["legacy_samples"].items() if k in sub["nodes"]}
        keep = {u["d0"] for u in used if u["d0"] in sub["nodes"]}
        cs2, _ = fit(sub, [u["file"] for u in used if u["d0"] in keep], mp_dps=120)
        flag = ("  [PRODUCTION contract]" if spec["production"] else "  [NOT a production contract]") + ("  LEGACY (frame unverified)" if any(u["legacy"] for u in used) else "")
        print(f"{name}: {len(used)} nodes{flag}")
        for k in range(-1, min(4, len(used) - 1)):
            c = cs[k]["c"]
            se = float(-log10(fabs(c - cs2[k]["c"]) / fabs(c))) if k in cs2 else None
            line = f"  c_{k:+d} = {nstr(c, 30)}" + (f"  (sub-grid {se:.1f} d)" if se else "")
            ex = exact[name].get(k)
            if ex is not None:
                line += f"  vs exact {float(-log10(fabs(c / ex() - 1))):.1f} d"
            elif k == 1:
                line += "  <- BLIND eps^1"
            line += f"  [+- {cs[k]['abs_err']}, amplification 10^{cs[k]['amplification_log10']:.2f}]{flag}"
            print(line)
        if "--out" in argv:
            out = argv[argv.index("--out") + 1] + f"-{name}.json"
            json.dump(seal({"schema": "pflow-fit.v1", "spec": spec, "used_samples": used,
                            "coefficients": {str(k): {x: c[x] for x in ("value", "abs_err", "amplification_log10")} for k, c in cs.items()}}),
                      open(out, "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        sys.exit(2)
