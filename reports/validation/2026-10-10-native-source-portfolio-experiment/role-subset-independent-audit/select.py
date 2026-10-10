#!/usr/bin/env python3
"""Read-only structural selection from the archived exact original rows."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
CORPUS = ROOT / "reports/validation/2026-10-09-finite-density-native-assembly/three-loop-singleton-source-diagnostic/original.json.gz"
GEOMETRY = ROOT / "reports/validation/2026-10-10-massless-reciprocal-energy-exploration/native-source-probe/geometry.json"
rows = json.loads(gzip.decompress(CORPUS.read_bytes()))["original_source_rows"]
roles = json.loads(GEOMETRY.read_text())["roles"]
arity = len(rows[0]["terms"][0]["label"])
roles += ["Ordinary"] * (arity - len(roles))
cut = [i for i, role in enumerate(roles) if role == "RequiredCut"]
occupation = [i for i, role in enumerate(roles) if role == "Occupation"]


def admitted(row, strict):
    # A literal index 0 is not the symbolic unchanged index a_i.
    return bool(row["terms"]) and all(
        all(term["label"][i] == [True, 0] for i in occupation)
        and all(
            term["label"][i][0]
            and (term["label"][i][1] == 0 if strict else term["label"][i][1] <= 0)
            for i in cut
        )
        for term in row["terms"]
    )


selections = {}
for name, strict in [("all_protected_shifts_zero", True), ("occupation_zero_cut_nonpositive", False)]:
    selected = [row for row in rows if admitted(row, strict)]
    selections[name] = {
        "count": len(selected),
        "original_ordinals": [row["original_ordinal"] for row in selected],
        "source_ids_for_inspection_only": [row["source_id"] for row in selected],
        "preserved_rows": selected,
    }
assert selections["all_protected_shifts_zero"]["count"] == 18
assert selections["occupation_zero_cut_nonpositive"]["count"] == 34
result = {
    "scope": "Read-only selector audit; no discovery, equations, source mutation or physical admission",
    "inputs_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (CORPUS, GEOMETRY)},
    "roles": roles,
    "original_source_count": len(rows),
    "selections": selections,
}
Path(__file__).with_name("selection.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({key: value["count"] for key, value in selections.items()}))
