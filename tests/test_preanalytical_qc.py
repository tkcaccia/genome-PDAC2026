import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "pipelines"
    / "preanalytical_qc"
    / "summarize_provider_nucleic_acid_qc.py"
)
SPEC = importlib.util.spec_from_file_location("provider_qc", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PreanalyticalQcTests(unittest.TestCase):
    def test_decimal_comma_and_missing_values(self):
        self.assertEqual(MODULE.number("42,3"), 42.3)
        self.assertIsNone(MODULE.number("---"))

    def test_ascat_group_summary(self):
        records = [
            {
                "ascat_boundary": True,
                "provider_qc": "PASS",
                "concentration_ng_ul": 50.0,
                "din": 3.0,
                "a260_280": 1.9,
                "a260_230": 2.1,
            },
            {
                "ascat_boundary": False,
                "provider_qc": "FAIL",
                "concentration_ng_ul": 25.0,
                "din": 2.5,
                "a260_280": 1.8,
                "a260_230": 2.3,
            },
        ]
        rows = MODULE.comparison_summary(records)
        pass_row = next(row for row in rows if row["metric"] == "provider_qc_pass")
        self.assertEqual(pass_row["boundary"], "1/1")
        self.assertEqual(pass_row["non_boundary"], "0/1")

    def test_coerce_qc_rows(self):
        rows = MODULE.coerce_qc_rows(
            [{"sample_id": "example", "din": "2,8", "provider_qc": " pass "}],
            ["din"],
        )
        self.assertEqual(rows[0]["din"], 2.8)
        self.assertEqual(rows[0]["provider_qc"], "PASS")


if __name__ == "__main__":
    unittest.main()
