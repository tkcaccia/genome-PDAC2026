#!/usr/bin/env python3
"""Site-aware adapters for isolated Nextflow execution of the validated scripts."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET


def run(command, root, env):
    print("Executing:", " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), cwd=root, env=env, check=True)


def environment(root):
    env = os.environ.copy()
    env["SLURM_SUBMIT_DIR"] = str(root)
    env["PYTHONPATH"] = str(root / "python_deps")
    env["MPLBACKEND"] = "Agg"
    env["PYTHON"] = sys.executable
    env["OMP_NUM_THREADS"] = os.environ.get("SLURM_CPUS_PER_TASK", "8")
    return env


def prepare(package, root, assets):
    package, root, assets = map(lambda value: Path(value).resolve(), (package, root, assets))
    if root == assets or root.exists():
        raise ValueError("Full runs require a new, nonexistent --run_root; existing results are protected")
    if root.name != "bed_2.5" or root.parent.name != "results":
        raise ValueError("--run_root must end in workspace/results/bed_2.5 for the legacy scripts")
    workspace = root.parents[1]
    if (workspace / "scripts").exists():
        raise ValueError("Use a dedicated empty workspace, not a pre-existing analysis workspace")
    root.mkdir(parents=True)
    shutil.copytree(package / "supporting_inputs/workspace", workspace, dirs_exist_ok=True)
    shutil.copytree(package / "scripts/legacy_sources", workspace / "scripts")
    shutil.copytree(package / "scripts", root / "scripts", ignore=shutil.ignore_patterns("legacy_sources", "__pycache__"))
    shutil.copytree(package / "raw_data/manifests", root / "manifests")
    shutil.copytree(package / "supporting_inputs/stale", root / "inputs/stale")
    shutil.copytree(package / "supporting_inputs/msisensor-pro-1.3.0", root / "inputs/msisensor-pro-1.3.0")
    for folder in ("logs", "alignments", "qc", "wes/calls", "wes/annotated", "svg_pdf"):
        (root / folder).mkdir(parents=True, exist_ok=True)
    with (package / "large_assets.csv").open() as handle:
        for row in csv.DictReader(handle):
            relative = Path(row["relative_path"])
            if relative.parts[0] == "alignments" or relative.suffix == ".list" or ".refseq." in relative.name:
                continue
            source = assets / relative
            if not source.is_file() or source.stat().st_size != int(row["bytes"]):
                raise ValueError(f"Missing or changed large input: {source}")
            destination = root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.symlink_to(source)
    deps = assets / "python_deps"
    if not deps.is_dir():
        raise FileNotFoundError(f"Recreate the recorded Python environment; missing {deps}")
    (root / "python_deps").symlink_to(deps)
    s10 = root / "svg_pdf/s10_rerun"
    s10.mkdir(parents=True)
    shutil.copy2(package / "svg_pdf/s10_rerun/generate_figure_s10_reconstruction.py", s10)
    tool = root / "inputs/msisensor-pro-1.3.0/binary/msisensor-pro-v1.3.0"
    tool.chmod(tool.stat().st_mode | 0o100)
    env = environment(root)
    run([sys.executable, root / "scripts/remap_bed_to_refseq.py",
         root / "inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.dict",
         root / "inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.bed",
         root / "inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.bed"], root, env)
    run(["java", "-jar", "/software/bio/picard/2.26.11/picard.jar", "BedToIntervalList",
         f"I={root}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.bed",
         f"O={root}/inputs/intervals/hg38_Twist_Bioscience_for_Illumina_Exome_2_5.refseq.interval_list",
         f"SD={root}/inputs/reference/GCF_000001405.40_GRCh38.p14_genomic.dict", "DROP_MISSING_CONTIGS=false"], root, env)


def audit(root, destination):
    root = Path(root).resolve()
    expected = {
        "bams": ("alignments", "*.dedup.bam", 28),
        "bam_indexes": ("alignments", "*.dedup.bam.bai", 28),
        "fastqc_html": ("qc/fastqc_2_5", "*_fastqc.html", 56),
        "coverage": ("qc/coverage_2_5", "*.hsmetrics.txt", 28),
        "corrected_hotspots": ("qc/coverage_2_5", "*.kras_hotspot_depth_v2.tsv", 28),
        "calls": ("wes/calls", "patient_*.somatic.vcf.gz", 14),
        "annotations": ("wes/annotated", "patient_*.vep.vcf.gz", 14),
        "annotation_indexes": ("wes/annotated", "patient_*.vep.vcf.gz.tbi", 14),
    }
    report = {"mode": "read-only completed-output audit", "root": str(root), "counts": {}}
    for label, (folder, pattern, count) in expected.items():
        files = list((root / folder).glob(pattern))
        if len(files) != count or any(path.stat().st_size == 0 for path in files):
            raise ValueError(f"{label}: expected {count} nonempty files, found {len(files)}")
        report["counts"][label] = len(files)
    with (root / "mofa/pdac_11/mofa/exports/mofa_factors.csv").open() as handle:
        factors = list(csv.DictReader(handle))
    if len(factors) != 11:
        raise ValueError("Primary MOFA model must contain 11 patients")
    with (root / "wes/msi_2_5/msi_scores_all_14.csv").open() as handle:
        msi = list(csv.DictReader(handle))
    if len(msi) != 14:
        raise ValueError("MSI table must contain all 14 patients")
    with (root / "svg_pdf/export_manifest.csv").open() as handle:
        exports = list(csv.DictReader(handle))
    for row in exports:
        for extension in ("svg", "pdf"):
            path = root / row[extension]
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(path)
            if extension == "svg":
                ET.parse(path)
    with (root / "svg_pdf/legacy_figure_coverage.csv").open() as handle:
        coverage = list(csv.DictReader(handle))
    if len(coverage) != 44:
        raise ValueError("Legacy coverage must account for all 44 old figures")
    report.update(primary_mofa_patients=11, msi_patients=14, exported_figure_pairs=len(exports),
                  legacy_figure_mappings=len(coverage), status="PASS",
                  scope="Artifact/cohort/export checks, not independent scientific validation")
    Path(destination).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def main():
    mode, *arguments = sys.argv[1:]
    if mode == "prepare":
        prepare(*arguments)
        return
    if mode == "audit":
        audit(*arguments)
        return
    root = Path(arguments[0]).resolve()
    env = environment(root)
    if mode == "batch":
        script, index = arguments[1:]
        env["SLURM_ARRAY_TASK_ID"] = index
        run(["bash", root / "scripts" / script], root, env)
    elif mode == "analyses":
        run(["bash", root / "scripts/21_full_analysis_2_5.sbatch"], root, env)
    elif mode == "msi_summary":
        run([sys.executable, root / "scripts/28_summarize_msi.py"], root, env)
        comparison = subprocess.run([sys.executable, str(root / "scripts/30_compare_msi_previous_2_5.py")], cwd=root, env=env)
        if comparison.returncode and not (root / "wes/msi_2_5/comparison_to_previous_msi_table.csv").is_file():
            comparison.check_returncode()
        if comparison.returncode:
            print("Old/new MSI differences recorded for scientific review; new scores are retained.")
    elif mode == "figures":
        run([sys.executable, root / "scripts/32_complete_stale_figure_inventory.py"], root, env)
        run([sys.executable, root / "scripts/37_figure7b_shared_pathways_2_5.py"], root, env)
        run([sys.executable, root / "scripts/36_figure7a_kras_status_2_5.py"], root, env)
    elif mode == "export":
        run([sys.executable, root / "svg_pdf/s10_rerun/generate_figure_s10_reconstruction.py"], root, env)
        run([sys.executable, root / "scripts/31_export_figures_svg_pdf_2_5.py"], root, env)
    else:
        raise ValueError(f"Unknown stage: {mode}")


if __name__ == "__main__":
    main()