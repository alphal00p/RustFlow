#!/usr/bin/env python3
"""Extract the effective Higgs-jet declarations without replacing the native SM.

This copies scientific model data, not model or tensor implementation. Run with
--check to verify that the two committed assets match their pinned source.
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "fixtures" / "gg-hg"
SOURCE_SHA256 = "1891f93a8b8926e1e66d54a95c03bba2d36e7ee1a812b91010fa0e88b4e92e12"


def encode(value):
    return (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode()


def extract():
    source_bytes = (DIRECTORY / "native-model.json").read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != SOURCE_SHA256:
        raise ValueError("The legacy scientific model does not match its pinned hash")
    source = json.loads(source_bytes)
    names = {
        "orders": {"GGGHEFT", "GGGHEW", "WW", "ZZ"},
        "parameters": {
            *(f"GGGH_HEFT_ForFac{i}" for i in range(1, 5)),
            *(f"GGGHEW{boson}_ForFac{i}_{part}"
              for boson in ("WW", "ZZ") for i in range(1, 5)
              for part in ("RE", "IM")),
            "CHEFT", "AllGGGHEWWWCoup", "AllGGGHEWZZCoup",
        },
        "couplings": {
            *(f"{prefix}_C{i}"
              for prefix in ("GGGH_HEFT", "GGGHEWWW", "GGGHEWZZ")
              for i in range(1, 5)),
        },
        "lorentz_structures": {f"gggHTensStruc{i}" for i in range(1, 5)},
        "vertex_rules": {"GGGHHEFT", "GGGHEWWW", "GGGHEWZZ"},
    }
    asset = {"schema": "rustflow-gg-hg-form-factor-vertices-v1"}
    for section, selected in names.items():
        records = [copy.deepcopy(item) for item in source[section]
                   if item["name"] in selected]
        if len(records) != len(selected) or {x["name"] for x in records} != selected:
            raise ValueError(f"Missing or duplicate scientific declarations: {section}")
        for item in records:
            if section == "parameters":
                item["value"] = [0.0, 0.0] if item["nature"] == "external" else None
                if item["name"] == "CHEFT":
                    if item["expression"] != "-1/3*aS/(𝜋*v)":
                        raise ValueError("Unexpected legacy HEFT normalization")
                    item["expression"] = "-1/3*aS/(𝜋*vev)"
            elif section == "couplings":
                item["value"] = None
            elif section == "vertex_rules" and item["particles"] != ["G", "G", "G", "H"]:
                raise ValueError("Unexpected legacy effective external states")
        asset[section] = records
    asset_bytes = encode(asset)
    provenance_bytes = (DIRECTORY / "native-model-provenance.json").read_bytes()
    source_provenance = json.loads(provenance_bytes)
    provenance = {
        "schema": "rustflow-gg-hg-form-factor-vertices-provenance-v1",
        "kind": "Effective interaction declarations extending the native Standard Model; no numerical integral boundaries or amplitudes",
        "source": {
            "path": "fixtures/gg-hg/native-model.json",
            "sha256": SOURCE_SHA256,
            "provenance_path": "fixtures/gg-hg/native-model-provenance.json",
            "provenance_sha256": hashlib.sha256(provenance_bytes).hexdigest(),
            "upstream": source_provenance["source"],
        },
        "asset": {
            "path": "fixtures/gg-hg/form-factor-vertices.json",
            "sha256": hashlib.sha256(asset_bytes).hexdigest(),
            "counts": {section: len(asset[section]) for section in names},
        },
        "transformations": [
            "Select only the three effective Higgs-jet vertices and their four orders, 23 parameters, 12 couplings and four Lorentz structures; preserve source order.",
            "Replace v by the existing Standard Model vev in CHEFT=-aS/(3*pi*vev); both native parameter definitions are exactly 2*MW*sw/ee.",
            "Set added internal-parameter and coupling cached values to null. Set external form-factor placeholders to zero; physical W/Z coefficients must be supplied and HEFT coefficients remain analytic kinematic functions.",
            "Keep legacy G/G/G/H labels in the declaration asset; the native extension must resolve their PDG21/25 species in the supplied model, without adding or renaming particles.",
        ],
        "tensor_ownership": "The four native UFO Lorentz expressions, color structures, vertex coupling slots and EW prefactors are copied unchanged. HEPKit/Spenso continue to generate diagrams and contract tensors; this asset implements no tensor algebra.",
        "normalization": "The integral projector retains -1/(mV2^2*(4*pi)^4) exactly once. These declarations add only the native coherent W/Z and HEFT coupling factors; no extra integral measure, spin/color average, flux or phase-space factor is introduced.",
        "regeneration": "python3 scripts/extract_gg_hg_form_factor_vertices.py; python3 scripts/extract_gg_hg_form_factor_vertices.py --check",
    }
    return {
        "form-factor-vertices.json": asset_bytes,
        "form-factor-vertices-provenance.json": encode(provenance),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Verify assets without writing them")
    args = parser.parse_args()
    for name, expected in extract().items():
        path = DIRECTORY / name
        if args.check:
            if not path.exists() or path.read_bytes() != expected:
                raise SystemExit(f"Stale or missing asset: {path}")
        else:
            path.write_bytes(expected)
        print(f"{'Checked' if args.check else 'Wrote'} {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
