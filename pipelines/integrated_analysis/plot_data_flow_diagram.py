#!/usr/bin/env python3
"""Draw a patient-data-safe multi-omics analysis workflow diagram."""

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
        "Multi-omics analysis workflow",
        ha="center",
        va="top",
        fontsize=19,
        fontweight="bold",
        family="serif",
    )
    ax.text(
        0.5,
        0.925,
        "Integration of matched tumour-normal WES and RNA-seq from African FFPE pancreatic and periampullary cancers",
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
        0.12,
        0.59,
        0.25,
        0.13,
        "Whole-exome sequencing",
        "nf-core/sarek 3.8.1\nMutect2 and Strelka2 SNV/indel calling\nASCAT CNA/purity; Manta SV",
        "#D9E6EC",
    )
    rna_box = add_box(
        ax,
        0.63,
        0.59,
        0.25,
        0.13,
        "RNA-seq processing",
        "nf-core/rnaseq 3.24.0\nSTAR counts and quality control\nGENCODE v46 gene annotation",
        "#DDE9DE",
    )
    add_arrow(ax, (0.44, 0.80), (0.245, 0.72))
    add_arrow(ax, (0.56, 0.80), (0.755, 0.72))

    variant_box = add_box(
        ax,
        0.02,
        0.36,
        0.21,
        0.13,
        "Somatic variants and drivers",
        "SNV/indel annotation\nKRAS hotspot classification\nTP53/CDKN2A/SMAD4 and DDR genes",
        "#D6E4EB",
    )
    msi_box = add_box(
        ax,
        0.265,
        0.36,
        0.21,
        0.13,
        "MSI and mutation burden",
        "MSIsensor-pro\nStrict rare-coding TMB\nMMR-gene evidence integration",
        "#E6DDD1",
    )
    de_box = add_box(
        ax,
        0.51,
        0.36,
        0.21,
        0.13,
        "Expression modelling",
        "Paired DESeq2 primary model\nedgeR quasi-likelihood\nlimma-voom sensitivity\nGene-level TPM derivation",
        "#D9E8D8",
    )
    tme_box = add_box(
        ax,
        0.755,
        0.36,
        0.21,
        0.13,
        "Pathway and TME scoring",
        "True GSVA and ssGSEA\nESTIMATE and MCP-counter\nEPIC, xCell and quanTIseq\nCohort-relative TME phenotypes",
        "#E1E6D5",
    )

    add_arrow(ax, (0.20, 0.59), (0.125, 0.49))
    add_arrow(ax, (0.29, 0.59), (0.37, 0.49))
    add_arrow(ax, (0.70, 0.59), (0.615, 0.49))
    add_arrow(ax, (0.80, 0.59), (0.86, 0.49))

    integrated_box = add_box(
        ax,
        0.21,
        0.13,
        0.58,
        0.12,
        "Integrated patient-level analysis",
        "Driver, MSI, TMB, CNA and SV features integrated with RNA programmes and immune/stromal scores\nExploratory genotype-phenotype associations and cohort-level interpretation",
        "#E4D9CD",
    )
    for x in (0.125, 0.37, 0.615, 0.86):
        add_arrow(ax, (x, 0.36), (0.50, 0.25))

    ax.text(
        0.5,
        0.055,
        "WES-derived CNA and SV inference has lower structural resolution than whole-genome sequencing; expression-defined phenotype contrasts are exploratory.",
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
