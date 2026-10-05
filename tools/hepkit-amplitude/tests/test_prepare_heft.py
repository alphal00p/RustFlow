import ast
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

script = Path(__file__).resolve().parents[1] / "scripts" / "prepare_heft.py"
spec = importlib.util.spec_from_file_location("prepare_heft", script)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ExactHeftInputs(unittest.TestCase):
    def test_closed_arithmetic_preserves_exact_inputs(self):
        expression = ast.parse("(pow(s, 2) - 3*t)/(2*s)", mode="eval").body
        self.assertEqual(module.evaluate(expression, {"s": Fraction(7, 3), "t": Fraction(-5, 11)}), Fraction(337, 231))

    def test_unsupported_or_float_expressions_are_rejected(self):
        for expression in ["0.5", "float(s)", "__import__('os')", "s.real", "pow(s, -1)", "pow(s, 99)"]:
            with self.subTest(expression=expression), self.assertRaises((ValueError, KeyError)):
                module.evaluate(ast.parse(expression, mode="eval").body, {"s": Fraction(2)})

    def fixture(self, base):
        # Synthetic scalar assignments exercise the generic reader; no upstream code.
        source = base / "bridge.cpp"
        source.write_text("\n".join(f"gghgHEFTTensor[{i}]={i+1}*s/t;" for i in range(4)))
        data = base / "input.json"
        data.write_text(json.dumps({"physical_coordinates": {"physical_s_t_MH_squared": ["7/3", "-5/11", "1"]}}))
        report = base / "report.json"
        report.write_text(json.dumps({"original": {"binary128": {"heft_squared": "1.234567890123456789e-3", "heft_ew_interference": "-2.234567890123456789e-5"}}}))
        output = base / "prepared"
        argv = [str(script), "--source", str(source), "--coherent-input", str(data), "--oracle-report", str(report), "--output", str(output)]
        return source, data, report, output, argv

    def test_preparation_preserves_strings_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            source, _, _, output, argv = self.fixture(Path(directory))
            with patch.object(module, "PINNED_BRIDGE_SHA256", hashlib.sha256(source.read_bytes()).hexdigest()), patch("sys.argv", argv):
                module.main()
                with self.assertRaisesRegex(ValueError, "overwrite"):
                    module.main()
            exact = json.loads((output / "exact-heft.json").read_text())
            self.assertEqual(exact["coefficients"][0]["value"], "-77/15")
            oracle = json.loads((output / "heft-oracle.json").read_text())
            self.assertEqual(oracle["binary128"]["heft_squared"], "1.234567890123456789e-3")

    def test_mismatched_external_source_creates_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, _, output, argv = self.fixture(Path(directory))
            with patch("sys.argv", argv), self.assertRaisesRegex(ValueError, "pinned"):
                module.main()
            self.assertFalse(output.exists())

    def test_binary_float_coordinate_creates_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            source, data, _, output, argv = self.fixture(Path(directory))
            data.write_text(json.dumps({"physical_coordinates": {"physical_s_t_MH_squared": [0.1, "-1", "1"]}}))
            with patch.object(module, "PINNED_BRIDGE_SHA256", hashlib.sha256(source.read_bytes()).hexdigest()), patch("sys.argv", argv), self.assertRaisesRegex(ValueError, "strings"):
                module.main()
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
