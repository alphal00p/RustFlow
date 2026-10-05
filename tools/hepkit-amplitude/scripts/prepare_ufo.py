#!/usr/bin/env python3
"""Prepare external pinned UFO input; never parse or contract tensor mathematics."""
from __future__ import annotations

import argparse
import ast
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import warnings

PINNED_REVISION = "89f64b93d0bdbd8ee90eb85737021189229b034a"


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False)
        stream.write("\n")


def read_json(path: Path) -> object:
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique)


def git(source: Path, *args: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(source), *args], stderr=subprocess.PIPE)


def snapshot(directory: Path) -> dict[str, str]:
    result = {}
    for path in sorted(directory.rglob("*")):
        if "__pycache__" in path.parts or path.suffix == ".pyc":
            continue
        if path.is_symlink():
            raise ValueError(f"symlinks are not supported in UFO input: {path}")
        if path.is_file():
            result[str(path.relative_to(directory))] = sha(path)
    return result


def calls(text: str, constructor: str) -> dict[str, dict[str, ast.expr]]:
    result = {}
    for node in ast.walk(ast.parse(text)):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Name):
            continue
        if node.func.id != constructor:
            continue
        keywords = {kw.arg: kw.value for kw in node.keywords}
        name = ast.literal_eval(keywords["name"])
        if not isinstance(name, str) or name in result:
            raise ValueError(f"invalid/duplicate {constructor} name: {name!r}")
        result[name] = keywords
    return result


def structures(text: str) -> dict[str, str]:
    result = {}
    for name, keywords in calls(text, "Lorentz").items():
        value = ast.literal_eval(keywords["structure"])
        if not isinstance(value, str):
            raise ValueError(f"nonliteral Lorentz structure: {name}")
        result[name] = value
    if not result:
        raise ValueError("no literal Lorentz definitions found")
    return result


def replace_nodes(text: str, replacements: list[tuple[ast.expr, str]]) -> str:
    # AST columns are UTF-8 byte offsets, including for non-ASCII source files.
    data = text.encode("utf-8")
    offsets = [0]
    for line in data.splitlines(keepends=True):
        offsets.append(offsets[-1] + len(line))
    edits = [(offsets[n.lineno - 1] + n.col_offset,
              offsets[n.end_lineno - 1] + n.end_col_offset, value.encode("utf-8"))
             for n, value in replacements]
    end = len(data)
    for start, stop, value in sorted(edits, reverse=True):
        if stop > end:
            raise ValueError("overlapping literal edits")
        data = data[:start] + value + data[stop:]
        end = start
    return data.decode("utf-8")


def prepare(args) -> None:
    source = args.source.resolve(strict=True)
    work = args.work.resolve()
    if work.exists() or work.is_relative_to(source):
        raise ValueError("work directory must be new and outside the source UFO")
    revision = git(source, "rev-parse", "HEAD").decode().strip()
    if revision != PINNED_REVISION:
        raise ValueError(f"expected pinned plugin {PINNED_REVISION}, got {revision}")
    repository = Path(git(source, "rev-parse", "--show-toplevel").decode().strip())
    relative = source.relative_to(repository)
    tracked = git(repository, "ls-tree", "-rz", "--name-only", "HEAD", "--", str(relative))
    expected = {}
    for raw in tracked.split(b"\0"):
        if not raw:
            continue
        name = raw.decode()
        expected[str(Path(name).relative_to(relative))] = hashlib.sha256(
            git(repository, "show", f"HEAD:{name}")).hexdigest()
    before = snapshot(source)
    if not expected or before != expected:
        raise ValueError("source UFO must exactly match its pinned tracked files (Python bytecode ignored)")
    texts = {name: (source / name).read_text(encoding="utf-8") for name in before if name.endswith(".py")}
    original_structures = structures(texts["lorentz.py"])
    changes = []

    def substitute(name: str, old: str, new: str, reason: str) -> None:
        if texts[name].count(old) != 1:
            raise ValueError(f"expected exactly one known compatibility occurrence in {name}: {old!r}")
        texts[name] = texts[name].replace(old, new)
        changes.append({"file": name, "old": old, "new": new, "reason": reason})

    substitute("object_library.py", ".iteritems()", ".items()", "Python 3 dictionary iteration")
    substitute("write_param_card.py", "print particles.__file__", "print(particles.__file__)", "Python 3 print syntax")
    substitute("write_param_card.py", "print 'write ./param_card.dat'", "print('write ./param_card.dat')", "Python 3 print syntax")
    substitute("__init__.py", "import particles\n", "import function_library\nimport particles\n", "Explicitly import registered function library before use")
    substitute("CT_couplings.py", "cmath.sqrt (2)", "cmath.sqrt(2)", "Function-call whitespace accepted by native UFO parser; identical Python expression")
    node = calls(texts["lorentz.py"], "Lorentz")["GHGHG"]["spins"]
    if ast.literal_eval(node) != [1, 1, 3]:
        raise ValueError("unexpected ghost spin signature")
    texts["lorentz.py"] = replace_nodes(texts["lorentz.py"], [(node, "[-1, -1, 3]")])
    changes.append({"file": "lorentz.py", "name": "GHGHG", "old": [1, 1, 3], "new": [-1, -1, 3], "reason": "Native HEPKit ghost spin signature; structure unchanged"})
    node = calls(texts["vertices.py"], "Vertex")["V_1001"]["lorentz"]
    if ast.dump(node) != ast.dump(ast.parse("[L.SSS1]", mode="eval").body):
        raise ValueError("unexpected four-scalar V_1001 metadata")
    texts["vertices.py"] = replace_nodes(texts["vertices.py"], [(node.elts[0], "L.SSSS1")])
    changes.append({"file": "vertices.py", "name": "V_1001", "old": "[L.SSS1]", "new": "[L.SSSS1]", "reason": "Four-scalar vertex requires four slots; both original scalar structures equal1"})
    lorentz = calls(texts["lorentz.py"], "Lorentz")
    if any(ast.literal_eval(lorentz[name]["structure"]) != "1" for name in ("SSS1", "SSSS1")):
        raise ValueError("scalar metadata normalization requires identical unit structures")
    if structures(texts["lorentz.py"]) != original_structures:
        raise ValueError("compatibility edits changed tensor expressions")
    with warnings.catch_warnings():
        # Existing UFO TeX labels contain backslash escapes; retain their bytes.
        warnings.simplefilter("ignore", SyntaxWarning)
        for name, text in texts.items():
            ast.parse(text, filename=name)
    work.mkdir(parents=True, exist_ok=False)
    destination = work / "UFO_model_gggH_python3"
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for name, text in texts.items():
        (destination / name).write_text(text, encoding="utf-8")
    if snapshot(source) != before:
        raise ValueError("source changed during preparation; discard this work directory")
    after = snapshot(destination)
    write_json(work / "original-lorentz.json", original_structures)
    write_json(work / "preparation.json", {
        "schema_version": 1, "source": str(source), "repository": str(repository), "source_revision": revision,
        "private_copy": str(destination), "python": sys.version, "script_sha256": sha(Path(__file__)),
        "changes": changes, "tensor_expressions_preserved": True,
        "files": {name: {"source": before[name], "prepared": after[name]} for name in before},
        "original_lorentz_sha256": sha(work / "original-lorentz.json"),
    })
    print(f"Prepared private Python 3 UFO and {len(original_structures)} unchanged Lorentz strings in {work}")


def apply_scopes(args) -> None:
    work = args.work.resolve(strict=True)
    provenance = read_json(work / "preparation.json")
    source = work / "UFO_model_gggH_python3"
    destination = work / "UFO_model_gggH_scoped"
    if destination.exists() or (work / "scoping-provenance.json").exists():
        raise ValueError("refusing to overwrite an existing scoped model/provenance")
    before = snapshot(source)
    if before != {name: hashes["prepared"] for name, hashes in provenance["files"].items()}:
        raise ValueError("prepared UFO files changed after recorded preparation")
    original_path = work / "original-lorentz.json"
    if sha(original_path) != provenance["original_lorentz_sha256"]:
        raise ValueError("original Lorentz input changed after preparation")
    original = read_json(original_path)
    changes = read_json(work / "scoped-lorentz.json")
    counts = read_json(work / "scope-counts.json")
    if not isinstance(changes, dict) or not changes or set(changes) != set(counts):
        raise ValueError("scope output requires matching nonempty structures and counts")
    if any(name not in original or not isinstance(value, str) or not value or
           not isinstance(counts[name], int) or isinstance(counts[name], bool) or counts[name] <= 0
           for name, value in changes.items()):
        raise ValueError("invalid or foreign scoped Lorentz names/strings/counts")
    text = (source / "lorentz.py").read_text(encoding="utf-8")
    if structures(text) != original:
        raise ValueError("prepared tensor strings differ from recorded raw input")
    nodes = calls(text, "Lorentz")
    text = replace_nodes(text, [(nodes[name]["structure"], repr(value)) for name, value in changes.items()])
    expected = dict(original, **changes)
    if structures(text) != expected:
        raise ValueError("scoped AST literal replacement mismatch")
    shutil.copytree(source, destination, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    (destination / "lorentz.py").write_text(text, encoding="utf-8")
    after = snapshot(destination)
    if [name for name in before if before[name] != after[name]] != ["lorentz.py"]:
        raise ValueError("unexpected non-Lorentz scoped-copy changes")
    write_json(work / "scoping-provenance.json", {
        "schema_version": 1, "source_revision": provenance["source_revision"], "private_source": str(source),
        "private_scoped_copy": str(destination), "preparation_sha256": sha(work / "preparation.json"),
        "script_sha256": sha(Path(__file__)), "original_lorentz_sha256": sha(original_path),
        "scoped_lorentz_sha256": sha(work / "scoped-lorentz.json"), "scope_counts_sha256": sha(work / "scope-counts.json"),
        "changed_names": sorted(changes), "counts": counts,
        "algorithm": "Native Symbolica Token/spenso::bracket scopes supplied by scope-lorentz; Python replaces only named AST structure string literals, with no tensor parsing or contraction.",
        "files": {name: {"before": before[name], "after": after[name]} for name in before},
    })
    print(f"Applied {len(changes)} named literal scopes to {destination}")


def extract_inputs(args) -> None:
    report_path = args.report.resolve(strict=True)
    output = args.output.resolve()
    if output.exists():
        raise ValueError("validation input directory must be new")
    report = read_json(report_path)
    scientific = report["scientific_input"]
    blocks = scientific["normalized_form_factors"]
    if {b["mass"] for b in blocks} != {"W", "Z"} or len(blocks) != 2:
        raise ValueError("expected one complete W and Z form-factor block")
    coordinates = scientific["physical_coordinates"]["physical_s_t_MH_squared"]
    if len(coordinates) != 3 or any(not isinstance(x, str) or not x for x in coordinates):
        raise ValueError("exact phase-point strings are required")
    factors = []
    for block in blocks:
        if sorted(f["index"] for f in block["form_factors"]) != [1, 2, 3, 4]:
            raise ValueError("expected four named form factors per mass")
        for value in block["form_factors"]:
            for number in (value["value"]["real"], value["value"]["imaginary"], value["absolute_error_allowance"]):
                if not isinstance(number, str) or not Decimal(number).is_finite():
                    raise ValueError("finite decimal strings required; binary64 values are not accepted")
            if Decimal(value["absolute_error_allowance"]) < 0:
                raise ValueError("negative form-factor error allowance")
        factors.append({"mass": block["mass"], "form_factors": [
            {"form_factor": f["index"], "normalized_value": f["value"],
             "normalized_absolute_error": f["absolute_error_allowance"]}
            for f in block["form_factors"]]})
    oracle = report["independent_oracle"]["result"]
    if not isinstance(oracle["binary128"]["SMATRIX"][0], str):
        raise ValueError("oracle decimal string missing")
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "amplitude-input.json", {"physical_coordinates": scientific["physical_coordinates"],
               "form_factors": factors, "normalization": scientific["normalization"]})
    write_json(output / "quad-result.json", oracle)
    write_json(output / "input-provenance.json", {"source_report": str(report_path), "source_report_sha256": sha(report_path),
               "mapping": "Exact JSON string/key mapping only; no numerical conversion, recomputation or precision promotion",
               "files": {name: sha(output / name) for name in ["amplitude-input.json", "quad-result.json"]}})
    print(f"Extracted recorded scientific inputs without float conversion to {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare", help="copy the pinned external UFO and record minimal compatibility edits")
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--work", type=Path, required=True)
    p.set_defaults(action=prepare)
    p = commands.add_parser("apply-scopes", help="apply native scope-lorentz output to a new private copy")
    p.add_argument("--work", type=Path, required=True)
    p.set_defaults(action=apply_scopes)
    p = commands.add_parser("extract-validation-inputs", help="extract exact input strings and original QP evidence from the durable ME report")
    p.add_argument("--report", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(action=extract_inputs)
    args = parser.parse_args()
    try:
        args.action(args)
    except (OSError, ValueError, KeyError, TypeError, subprocess.CalledProcessError, SyntaxError, InvalidOperation, IndexError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
