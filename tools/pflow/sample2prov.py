"""Add rows to a provided table (pflow-provided.v2) straight from ACCEPTED route-B samples of the same d0: the massless value of an
integer-index integral of a family graph, E-normalised, as the sample gives it (a master at m = 1, or an --extra target), with the
sample's achieved abs_err. Used for lower-loop integrals that the p-integral answers do not contain (or as an independent source).
usage: sample2prov.py --table IN.json --out OUT.json --add GRAPH.json INDEX SAMPLE.json FRAME_SHA256 [--add ...] [--measurement]
  GRAPH.json: the family graph (lines[].endpoints, external.in/out); INDEX: comma list of the line indices; FRAME_SHA256: the
  sample's required frame digest.
Refused (exit 2): an unsealed or edited table, a table or sample of another d0, a rejected, verdict-less or edited sample, a frame
digest mismatch, a downgraded sample (no two-precision pass, unreachable inputs) unless --measurement (then recorded in the table's
downgrades), a target the sample does not contain, a row that duplicates an existing one beyond the two uncertainties. The output
is sealed again."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pflowio import Refused, PROVIDED_SCHEMA, rat, read_sample, load_json, sha256_file, rel_digits, check_seal, seal


def main(argv):
    from mpmath import mp, mpf, fabs
    mp.dps = 600
    opt = lambda k: argv[argv.index(k) + 1] if k in argv else None
    t = load_json(opt("--table"))
    if t.get("schema") != PROVIDED_SCHEMA:
        raise Refused(f"{opt('--table')}: schema {t.get('schema')!r}")
    check_seal(opt("--table"), t)
    d0 = rat(t["d0"])
    adds = [argv[i + 1:i + 5] for i, a in enumerate(argv) if a == "--add"]
    if any(len(a) != 4 for a in adds):
        raise Refused("--add GRAPH.json INDEX SAMPLE.json FRAME_SHA256")
    if not adds or not opt("--out"):
        raise Refused(__doc__.split("usage: ")[1].split("\n")[0])
    for gpath, idxs, spath, fsha in adds:
        g = load_json(gpath)
        idx = [int(x) for x in idxs.split(",")]
        lines = g["lines"]
        # the family index may carry irreducible numerators after the lines: only absent ones (0) are allowed here
        if len(idx) < len(lines) or any(x != 0 for x in idx[len(lines):]):
            raise Refused(f"{gpath}: {len(lines)} lines, index {idx} (numerators are not supported)")
        s = read_sample(spath, d0=d0, frame_sha=fsha)
        if s["downgrades"]:
            if "--measurement" not in argv:
                raise Refused(f"{spath}: downgraded sample ({'; '.join(s['downgrades'])}); --measurement records it instead")
            t.setdefault("downgrades", []).extend(f"row {idxs}: {x}" for x in s["downgrades"])
        tgt = tuple(idx)
        if tgt in s["extras"]:
            v, e = s["extras"][tgt]
        elif (1, tgt) in s["values"]:
            v, e = s["values"][(1, tgt)]
        else:
            raise Refused(f"{spath}: no value of {idx}")
        edges = [[l["endpoints"][0], l["endpoints"][1], a] for l, a in zip(lines, idx[:len(lines)])]
        fam = Path(gpath).stem.removeprefix("graph-")
        row = {"edges": edges, "ext": [g["external"]["in"], g["external"]["out"]], "value": v, "abs_err": e,
               "id": f"sample:{fam}:{idxs}", "family": fam, "index": idx, "source": Path(spath).name}
        for r in t["values"]:
            if r.get("family") == row["family"] and r.get("index") == idx:
                if fabs(mpf(r["value"]) - mpf(v)) > mpf(r["abs_err"]) + mpf(e):
                    raise Refused(f"{spath}: {idx} disagrees with the table row {r.get('id')} beyond the uncertainties")
        t["values"].append(row)
        t["sources"].append({"file": str(spath), "sha256": s["sha256"], "role": f"sample for {row['family']} {idx}"})
        dg = rel_digits(v, e)
        old = float(t["worst_relative_digits"]) if t.get("worst_relative_digits") is not None else dg
        t["worst_relative_digits"] = f"{min(old, dg):.2f}"
    t["producer"] = t.get("producer", "") + " + RustFlow tools/pflow/sample2prov.py"
    seal(t)
    with open(opt("--out"), "w") as f:
        json.dump(t, f, indent=1, ensure_ascii=False)
    print(len(adds), "rows added; table", len(t["values"]), "rows, worst", t["worst_relative_digits"], "relative digits")


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except Refused as x:
        print(f"REFUSED: {x}", file=sys.stderr)
        sys.exit(2)
