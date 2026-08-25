#!/usr/bin/env python3
"""Patient-data-free regression tests for the public workflow interfaces."""

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


class PublicWorkflowTests(unittest.TestCase):
    def run_script(self, script: Path, *arguments: object) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["python3", str(script), *(str(value) for value in arguments)],
            cwd=ROOT,
            check=True,
            text=True,
            capture_output=True,
        )

    def test_example_configuration_renders(self) -> None:
        result = self.run_script(
            ROOT / "scripts" / "render_config_env.py",
            ROOT / "config" / "config.example.yaml",
        )
        self.assertIn("RUN_STANDARD_DE=true", result.stdout)
        self.assertIn("RUN_PHENOTYPE_ASSIGNMENT=false", result.stdout)

    def test_score_assembly_and_ranked_phenotype_assignment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            sample_ids = [f"EXAMPLE_{index:02d}" for index in range(1, 15)]

            immune = pd.DataFrame(
                {
                    "feature": ["immune score", "T cell CD8+"],
                    **{
                        sample: [8.0, 7.0] if index <= 3 else ([1.0, 1.0] if index <= 6 else [4.0, 4.0])
                        for index, sample in enumerate(sample_ids, start=1)
                    },
                }
            )
            programmes = pd.DataFrame(
                {
                    "programme": ["CAF_ECM", "EMT_INVASION"],
                    **{
                        sample: [1.0, 1.0] if index <= 3 else ([8.0, 7.0] if index <= 6 else [4.0, 4.0])
                        for index, sample in enumerate(sample_ids, start=1)
                    },
                }
            )
            metadata = pd.DataFrame(
                {
                    "sample_id": sample_ids,
                    "patient_id": [f"EXAMPLE_PATIENT_{index:02d}" for index in range(1, 15)],
                    "condition": ["Tumour"] * 14,
                }
            )

            immune_file = temporary / "immune.tsv"
            programme_file = temporary / "programmes.tsv"
            metadata_file = temporary / "metadata.tsv"
            assembled_file = temporary / "assembled.tsv"
            prefix = temporary / "phenotype"
            immune.to_csv(immune_file, sep="\t", index=False)
            programmes.to_csv(programme_file, sep="\t", index=False)
            metadata.to_csv(metadata_file, sep="\t", index=False)

            self.run_script(
                ROOT / "pipelines" / "phenotype_assignment" / "assemble_tme_score_table.py",
                "--score-table",
                f"estimate={immune_file}",
                "--programme-table",
                programme_file,
                "--metadata",
                metadata_file,
                "--output",
                assembled_file,
            )
            assembled = pd.read_csv(assembled_file, sep="\t")
            self.assertEqual(len(assembled), 14)
            self.assertIn("estimate_immune_score", assembled.columns)
            self.assertIn("programme_emt_invasion", assembled.columns)

            self.run_script(
                ROOT / "pipelines" / "phenotype_assignment" / "assign_tme_phenotype_groups.py",
                "--scores",
                assembled_file,
                "--immune-columns",
                "estimate_immune_score,estimate_t_cell_cd8_plus",
                "--stromal-columns",
                "programme_caf_ecm",
                "--emt-columns",
                "programme_emt_invasion",
                "--target-per-extreme",
                3,
                "--out-prefix",
                prefix,
            )
            assignments = pd.read_csv(prefix.with_suffix(".assignments.tsv"), sep="\t")
            counts = assignments["phenotype_group"].value_counts().to_dict()
            self.assertEqual(counts["ImmuneHigh_StromalLow"], 3)
            self.assertEqual(counts["StromalHigh_EMTHigh_ImmuneLow"], 3)
            self.assertEqual(counts["Intermediate_or_mixed"], 8)
            self.assertEqual(assignments.loc[:2, "phenotype_group"].nunique(), 1)
            self.assertEqual(assignments.loc[3:5, "phenotype_group"].nunique(), 1)


if __name__ == "__main__":
    unittest.main()
