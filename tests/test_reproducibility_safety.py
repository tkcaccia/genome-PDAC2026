"""Synthetic-only tests: no patient files, fitted models or scientific tools."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from scripts.mofa.model_inputs import align_matrix, require_axis, write_input_identity

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("workflow_stage", ROOT / "workflow/nextflow/stage.py")
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)

try:
    import pandas as pd
except ImportError:
    pd = None


class AxisTests(unittest.TestCase):
    def test_saved_order_not_observed_order(self):
        self.assertEqual(require_axis(["beta", "alpha"], ["alpha", "beta"], "samples"),
                         ["alpha", "beta"])

    def test_missing_and_duplicate_names_fail(self):
        with self.assertRaises(ValueError):
            require_axis(["alpha", "alpha"], ["alpha"], "samples")
        with self.assertRaises(ValueError):
            require_axis(["alpha"], ["alpha", "beta"], "samples")
        with self.assertRaises(ValueError):
            require_axis(["alpha"], [], "samples")

    def test_extra_samples_fail_but_unselected_qc_features_are_allowed(self):
        with self.assertRaises(ValueError):
            require_axis(["alpha", "beta"], ["alpha"], "samples")
        self.assertEqual(require_axis(["gene_a", "gene_b"], ["gene_a"], "features", True),
                         ["gene_a"])


@unittest.skipIf(pd is None, "Install the test-only numpy/pandas dependencies")
class MatrixTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame([[1, 0], [0, 1], [1, 1]],
                            index=["gene_b", "gene_a", "unselected"],
                            columns=["sample_b", "sample_a"])

    def test_reorders_without_changing_or_restandardising_values(self):
        out = align_matrix(self.frame(), ["gene_a", "gene_b"],
                           ["sample_a", "sample_b"], "wes")
        self.assertEqual(out.values.tolist(), [[1.0, 0.0], [0.0, 1.0]])

    def test_rejects_missing_nonfinite_and_nonbinary_values(self):
        frame = self.frame()
        with self.assertRaises(ValueError):
            align_matrix(frame, ["missing"], list(frame.columns), "wes")
        frame.iloc[0, 0] = float("nan")
        with self.assertRaises(ValueError):
            align_matrix(frame, ["gene_b"], list(frame.columns), "wes")
        frame.iloc[0, 0] = 0.5
        with self.assertRaises(ValueError):
            align_matrix(frame, ["gene_b"], list(frame.columns), "wes")

    def test_identity_does_not_overwrite_existing_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "identity.json"
            frame = self.frame()
            reference = pd.DataFrame([[0.1], [0.2]], index=list(frame.columns),
                                     columns=["factor_1"])
            write_input_identity(target, {"wes": frame}, reference)
            before = target.read_bytes()
            with self.assertRaises(FileExistsError):
                write_input_identity(target, {"wes": frame}, reference)
            self.assertEqual(target.read_bytes(), before)


class WorkflowSafetyTests(unittest.TestCase):
    def test_public_repository_is_not_a_private_destination(self):
        with self.assertRaises(ValueError):
            stage.private_location(ROOT / "new_output", "Run")

    def test_rejects_path_traversal_and_absolute_assets(self):
        for value in ("../escape", "/absolute/input", "nested/../../escape", "."):
            with self.assertRaises(ValueError):
                stage.relative_asset(value)
        self.assertEqual(stage.relative_asset("inputs/reference.fna"),
                         Path("inputs/reference.fna"))

    def test_existing_workspace_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            workspace = base / "workspace"
            workspace.mkdir()
            important = workspace / "keep.txt"
            important.write_text("retain")
            with self.assertRaises(ValueError):
                stage.validate_roots(base / "package", workspace / "results/bed_2.5",
                                     base / "assets")
            self.assertEqual(important.read_text(), "retain")

    def test_inputs_and_outputs_cannot_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            with self.assertRaises(ValueError):
                stage.validate_roots(base / "package", base / "assets/results/bed_2.5",
                                     base / "assets")

    def test_preparation_uses_only_synthetic_assets_and_mocked_commands(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            package, assets = base / "package", base / "assets"
            for relative in ("supporting_inputs/workspace", "scripts/legacy_sources",
                             "raw_data/manifests", "supporting_inputs/stale",
                             "supporting_inputs/msisensor-pro-1.3.0/binary", "svg_pdf/s10_rerun"):
                (package / relative).mkdir(parents=True, exist_ok=True)
            (package / "large_assets.csv").write_text("relative_path,bytes\n")
            (package / "supporting_inputs/msisensor-pro-1.3.0/binary/msisensor-pro-v1.3.0").write_text("synthetic")
            (package / "svg_pdf/s10_rerun/generate_figure_s10_reconstruction.py").write_text("# synthetic\n")
            (assets / "python_deps").mkdir(parents=True)
            destination = base / "new_workspace/results/bed_2.5"
            with patch.object(stage, "run") as mocked:
                stage.prepare(package, destination, assets)
            self.assertEqual(mocked.call_count, 2)
            self.assertEqual((destination / "scripts/26_msi_scan.sbatch").read_bytes(),
                             (ROOT / "scripts/wes/raw/26_msi_scan.sbatch").read_bytes())

    def test_audit_cannot_write_inside_its_input_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError):
                stage.audit(root, root / "report.json")

    def test_missing_asset_fails_before_creating_workspace(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            package = base / "package"
            package.mkdir()
            (package / "large_assets.csv").write_text("relative_path,bytes\ninputs/missing.bin,99\n")
            destination = base / "new_workspace/results/bed_2.5"
            with self.assertRaises(ValueError):
                stage.prepare(package, destination, base / "assets")
            self.assertFalse((base / "new_workspace").exists())


if __name__ == "__main__":
    unittest.main()
