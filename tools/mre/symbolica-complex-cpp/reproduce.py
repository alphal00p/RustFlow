"""Compile and evaluate i*x through the public scalar Symbolica Python API."""

import argparse
from contextlib import nullcontext
import json
from pathlib import Path
import tempfile

from symbolica import E, S


def check_values(values, label):
    expected = [-3 + 2j, -4j, 5 + 0j]
    actual = [complex(row[0]) for row in values]
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, received {actual!r}")
    return [str(value) for value in actual]


def reproduce(directory):
    directory.mkdir(parents=True, exist_ok=True)
    x = S("complex_cpp_mre::x")
    inputs = [[2 + 3j], [-4 + 0j], [-5j]]
    results = []
    for direct in (True, False):
        name = "complex_direct" if direct else "complex_linear"
        record = {"direct_translation": direct}
        evaluator = (E("1i") * x).evaluator(
            [x], jit_compile=False, direct_translation=direct
        )
        # This independent interpreter control must already produce i*x.
        record["interpreted_values"] = check_values(
            evaluator.evaluate_complex(inputs), name + " interpreted"
        )
        cpp = directory / (name + ".cpp")
        try:
            compiled = evaluator.compile(
                name,
                str(cpp),
                str(directory / (name + ".so")),
                number_type="complex",
                inline_asm="none",
                native=False,
            )
            record["compiled_values"] = check_values(
                compiled.evaluate(inputs), name + " compiled"
            )
        except (ValueError, AssertionError) as error:
            record["status"] = "failed"
            record["error"] = str(error)
            (directory / (name + ".error")).write_text(str(error))
        else:
            record["status"] = "passed"
        record["generated_assignments"] = [
            line.strip()
            for line in cpp.read_text().splitlines()
            if "params[0]" in line
        ]
        results.append(record)
    (directory / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    return 1 if any(record["status"] == "failed" for record in results) else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-directory",
        type=Path,
        help="Keep generated C++, libraries, and results here; otherwise use a temporary directory.",
    )
    args = parser.parse_args()
    context = (
        nullcontext(args.output_directory)
        if args.output_directory is not None
        else tempfile.TemporaryDirectory(prefix="symbolica-complex-cpp-")
    )
    with context as directory:
        return reproduce(Path(directory).resolve())


if __name__ == "__main__":
    raise SystemExit(main())
