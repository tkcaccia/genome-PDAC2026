import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "pipelines" / "sarek_sv_cna" / "audit_ascat_purity.py"


class AscatPurityAuditTests(unittest.TestCase):
    def test_boundary_fit_is_flagged_and_summarized(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            molecular = tmpdir / "molecular.tsv"
            clinical = tmpdir / "clinical.tsv"
            output = tmpdir / "output"

            molecular_rows = [
                {
                    "patient_id": "A",
                    "ASCAT_aberrant_cell_fraction_purity_estimate": "1.000",
                    "ASCAT_ploidy": "2.001",
                    "ASCAT_purity_ploidy_source": "purityploidy.txt",
                    "large_CNV_event_count": "1",
                    "cnv_burden_altered_segments": "1",
                    "segments": "24",
                },
                {
                    "patient_id": "B",
                    "ASCAT_aberrant_cell_fraction_purity_estimate": "0.45",
                    "ASCAT_ploidy": "2.1",
                    "ASCAT_purity_ploidy_source": "purityploidy.txt",
                    "large_CNV_event_count": "8",
                    "cnv_burden_altered_segments": "8",
                    "segments": "31",
                },
                {
                    "patient_id": "C",
                    "ASCAT_aberrant_cell_fraction_purity_estimate": "0.70",
                    "ASCAT_ploidy": "2.8",
                    "ASCAT_purity_ploidy_source": "purityploidy.txt",
                    "large_CNV_event_count": "20",
                    "cnv_burden_altered_segments": "20",
                    "segments": "45",
                },
            ]
            with molecular.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(molecular_rows[0]), delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(molecular_rows)

            clinical_rows = [
                {"patient_id": "A", "Tumor_stroma_ratio": "40/60"},
                {"patient_id": "B", "Tumor_stroma_ratio": "50/50"},
                {"patient_id": "C", "Tumor_stroma_ratio": "70/30"},
            ]
            with clinical.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=list(clinical_rows[0]), delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(clinical_rows)

            subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--molecular-table",
                    str(molecular),
                    "--clinical-table",
                    str(clinical),
                    "--output-dir",
                    str(output),
                ],
                check=True,
                capture_output=True,
                text=True,
            )

            summary = (output / "ascat_purity_audit_summary.md").read_text()
            self.assertIn("Exact upper-boundary estimates (purity >= 0.9995): 1", summary)
            self.assertIn("Boundary estimates with near-diploid ploidy (1.90-2.10): 1/1", summary)

            with (output / "ascat_purity_patient_level_audit.tsv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["pathology_tumour_fraction"], "0.4")
            self.assertIn("purity_at_upper_boundary", rows[0]["diagnostic_flags"])
            self.assertIn("discordant_with_pathology_fraction", rows[0]["diagnostic_flags"])


if __name__ == "__main__":
    unittest.main()
