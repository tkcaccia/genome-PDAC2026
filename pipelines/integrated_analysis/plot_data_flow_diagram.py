#!/usr/bin/env python3
"""Draw a patient-data-safe analysis and provenance diagram."""

from __future__ import annotations

import argparse
from pathlib import Path

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort-size", type=int, default=None)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def load_plotting_dependencies() -> None:
    global plt, FancyArrowPatch, FancyBboxPatch
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required to draw the data-flow figure. Install the "
            "environment documented in env/environment.yml."
        ) from exc


def add_box(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    colour: str,
    linestyle: str = "solid",
):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.5,
        edgecolor="#363636",
        facecolor=colour,
        linestyle=linestyle,
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height * 0.68,
        title,
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        family="serif",
    )
    ax.text(
        x + width / 2,
        y + height * 0.34,
        body,
        ha="center",
        va="center",
        fontsize=8.6,
        family="serif",
        linespacing=1.25,
    )
    return patch


def add_arrow(ax, start, end, dashed: bool = False):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=13,
        linewidth=1.25,
        color="#555555",
        linestyle="dashed" if dashed else "solid",
        connectionstyle="arc3,rad=0",
    )
    ax.add_patch(arrow)


def main() -> None:
    args = parse_args()
    load_plotting_dependencies()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cohort_text = (
        f"{args.cohort_size} FFPE pancreatic/periampullary cases"
        if args.cohort_size is not None
        else "FFPE pancreatic/periampullary cohort"
    )

    fig, ax = plt.subplots(figsize=(15.5, 9.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(
        0.5,
        0.965,
        "Multi-omics analysis flow and current provenance status",
        ha="center",
        va="top",
        fontsize=19,
        fontweight="bold",
        family="serif",
    )
    ax.text(
        0.5,
        0.925,
        "Solid borders indicate outputs available in the canonical audit tree; dashed borders indicate historical records without recovered primary outputs",
        ha="center",
        va="top",
        fontsize=9.5,
        family="serif",
        color="#444444",
    )

    input_box = add_box(
        ax,
        0.35,
        0.80,
        0.30,
        0.09,
        "FFPE tumour and matched normal material",
        f"{cohort_text}\nWES and paired RNA-seq where available",
        "#E8E1D2",
    )

    wes_box = add_box(
        ax,
        0.06,
        0.59,
        0.25,
        0.13,
        "WES processing and summaries",
        "nf-core/sarek 3.8.1 record\nMutect2/Strelka SNV-indel review\nASCAT CNA/purity; Manta SV",
        "#D9E6EC",
    )
    rna_box = add_box(
        ax,
        0.375,
        0.59,
        0.25,
        0.13,
        "RNA-seq processing",
        "nf-core/rnaseq 3.24.0\nSTAR counts and quality control\nGENCODE v46 gene annotation",
        "#DDE9DE",
    )
    fusion_box = add_box(
        ax,
        0.69,
        0.59,
        0.25,
        0.13,
        "RNA-fusion workflow record",
        "nf-core/rnafusion 4.1.0 historical record\nPrimary run logs/results not recovered\nNot used for biological claims",
        "#F1E6D4",
        linestyle="dashed",
    )

    add_arrow(ax, (0.44, 0.80), (0.185, 0.72))
    add_arrow(ax, (0.50, 0.80), (0.50, 0.72))
    add_arrow(ax, (0.56, 0.80), (0.815, 0.72), dashed=True)

    variant_box = add_box(
        ax,
        0.03,
        0.36,
        0.21,
        0.13,
        "Genomic review",
        "COSMIC annotation summaries\nKRAS coverage/variant audit\nTP53/CDKN2A/SMAD4 integration",
        "#D6E4EB",
    )
    msi_box = add_box(
        ax,
        0.27,
        0.36,
        0.21,
        0.13,
        "MSI and mutation burden",
        "MSIsensor-pro across all pairs\nStrict rare-coding TMB\nCaller-consensus sensitivity\nMMR-gene evidence integration",
        "#E6DDD1",
    )
    de_box = add_box(
        ax,
        0.52,
        0.36,
        0.21,
        0.13,
        "Expression modelling",
        "Paired DESeq2 primary model\nedgeR quasi-likelihood\nlimma-voom sensitivity\nGene-level TPM derivation",
        "#D9E8D8",
    )
    tme_box = add_box(
        ax,
        0.76,
        0.36,
        0.21,
        0.13,
        "Pathway and TME scoring",
        "True GSVA and ssGSEA\nESTIMATE and MCP-counter\nEPIC, xCell and quanTIseq\nCohort-relative TME phenotypes",
        "#E1E6D5",
    )

    add_arrow(ax, (0.14, 0.59), (0.135, 0.49))
    add_arrow(ax, (0.24, 0.59), (0.375, 0.49))
    add_arrow(ax, (0.44, 0.59), (0.625, 0.49))
    add_arrow(ax, (0.56, 0.59), (0.865, 0.49))

    integrated_box = add_box(
        ax,
        0.21,
        0.13,
        0.58,
        0.12,
        "Integrated, patient-level interpretation",
        "Anonymized driver/MSI/TMB/CNA/SV summaries + RNA programmes + immune/stromal features\nExploratory associations, manuscript figures, tables and auditable method notes",
        "#E4D9CD",
    )
    for x in (0.135, 0.375, 0.625, 0.865):
        add_arrow(ax, (x, 0.36), (0.50, 0.25))

    ax.text(
        0.5,
        0.055,
        "Interpretive boundary: WES-compatible CNA/SV inference is lower resolution than WGS; expression-defined phenotype contrasts are descriptive; no private data are published in the code repository.",
        ha="center",
        va="center",
        fontsize=9,
        family="serif",
        color="#474747",
    )

    fig.savefig(args.outdir / "multiomics_data_flow_provenance.png", dpi=320, bbox_inches="tight")
    fig.savefig(args.outdir / "multiomics_data_flow_provenance.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
