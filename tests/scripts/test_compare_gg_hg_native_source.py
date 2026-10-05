"""Comparator mechanics only: synthetic reports never enter a numerical bank."""

from copy import deepcopy
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "compare_gg_hg_native_source", ROOT / "scripts/compare_gg_hg_native_source.py"
)
COMPARISON = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARISON)


class NativeSourceComparison(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.native_path = Path(self.directory.name) / "synthetic-native.json"
        self.reference_path = Path(self.directory.name) / "synthetic-reference.json"
        values = [[{"real": "0", "imaginary": "0"} for _ in range(48)] for _ in range(5)]
        values[4][17] = {"real": "1.2345678901234567890123456789", "imaginary": "-2.345678901234567890123456789"}
        self.reference = {
            "provenance": {"scope": "synthetic comparator unit test, not scientific evidence"},
            "systems": [{"family": "Planar_EW1", "dimension": 48, "coordinate_names": ["s", "t", "b"], "roots": [{"name": "root1"}, {"name": "root2"}]}],
            "cases": [{"label": "W-Planar_EW1-1", "family": "Planar_EW1", "dimension": 48,
                       "start": ["10", "-3", "2"], "source_root_germs": ["principal", "opposite_principal"],
                       "source_values": values, "source_errors": [["1e-26"] * 48 for _ in range(5)],
                       "source_verified_digits_cap": 24, "source_sha256": "synthetic",
                       "source_delta": "1e-27", "source_error_decimal_exponent": -26, "matrix_sha256": "synthetic"}],
        }
        self.native = {
            "schema": "rustflow-native-gg-hg-boundary-run-v1", "status": "success", "verified_bank_saved": True,
            "family": "Planar_EW1", "label": "W-Planar_EW1-1", "namespace": "test",
            "mathematical_source": "https://arxiv.org/abs/2112.07578v1",
            "physical_map_source_sha256": COMPARISON.PHYSICAL_MAP_SOURCES["Planar_EW1"],
            "canonical_identity": "synthetic", "build": {"scope": "unit test only"},
            "coordinates": {"test::{}::s": "10", "test::{}::t": "-3", "test::{}::b": "2"},
            "root_germs": {"test::{}::root1": "Principal", "test::{}::root2": "Opposite"},
            "boundary": {"identity": "synthetic", "provenance": "Native automatic auxiliary-mass boundary; extracted-map synthetic; no numerical reference seed; SYNTHETIC TEST",
                         "leading_epsilon_power": 0, "last_epsilon_power": 4, "verified_digits": 30,
                         "canonical_basis_indices": list(range(1, 49)), "coefficients": deepcopy(values),
                         "absolute_errors": [["1e-35"] * 48 for _ in range(5)]},
        }

    def compare(self, digits=20):
        self.native_path.write_text(json.dumps(self.native))
        self.reference_path.write_text(json.dumps(self.reference))
        return COMPARISON.compare_files(self.native_path, self.reference_path, digits)

    def test_complete_match_preserves_reference_cap(self):
        result = self.compare()
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["checked_coefficients"], 240)
        self.assertEqual(result["reference_verified_digits_cap"], 24)
        self.assertEqual(result["native_verified_digits"], 30)
        self.assertEqual(result["corroborated_mixed_digits"], 20)
        self.assertTrue(result["reference_opened_after_native_success"])

    def test_failed_native_run_never_opens_reference(self):
        self.native["status"] = "prepared"
        self.native_path.write_text(json.dumps(self.native))
        with self.assertRaisesRegex(ValueError, "references were not opened"):
            COMPARISON.compare_files(self.native_path, self.reference_path)
        self.assertFalse(self.reference_path.exists())

    def test_recorded_error_is_stronger_than_requested_twenty_digits(self):
        self.native["boundary"]["coefficients"][2][3]["imaginary"] = "1.000000002e-26"
        result = self.compare()
        self.assertEqual(result["status"], "reference_disagreement")
        failure = result["failures"][0]
        self.assertEqual((failure["epsilon_power"], failure["canonical_basis_index"]), (2, 4))
        self.assertFalse(failure["within_combined_recorded_errors"])
        self.assertTrue(failure["within_requested_mixed_accuracy"])

    def test_exact_zero_has_no_invented_reference_error(self):
        self.reference["cases"][0]["source_errors"][0][0] = "0"
        self.native["boundary"]["absolute_errors"][0][0] = "0"
        self.native["boundary"]["coefficients"][0][0]["real"] = "1e-120"
        result = self.compare()
        self.assertEqual(result["status"], "reference_disagreement")
        self.assertEqual(result["maximum_ratio_to_combined_recorded_errors"], "infinity")

    def test_mismatched_point_and_sheet_are_rejected(self):
        self.native["coordinates"]["test::{}::s"] = "10001/1000"
        with self.assertRaisesRegex(ValueError, "exact reference SOURCE"):
            self.compare()
        self.native["coordinates"]["test::{}::s"] = "10"
        self.native["root_germs"]["test::{}::root2"] = "Principal"
        with self.assertRaisesRegex(ValueError, "root germs"):
            self.compare()

    def test_mantissa_cannot_upgrade_precision(self):
        with self.assertRaisesRegex(ValueError, "cap 24"):
            self.compare(digits=25)


if __name__ == "__main__":
    unittest.main()
