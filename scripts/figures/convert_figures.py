#!/usr/bin/env python3
# Current collaborator code export. Final supplied vector presentation scripts, consuming private derived tables.
# Source: PDAC_Manuscript_Figs_Sne/convert_figures.py
# Private paths/cohort labels are supplied separately; see docs/inputs.md.
# Publication export only: not run or scientifically revalidated during reorganization.
"""Create SVG counterparts and preserve PDF outputs for supplied figure PDFs."""

from __future__ import annotations

from pdac_runtime import project_root, figure_root, private_path, private_value, cohort as configured_cohort, code_path

import shutil
import argparse
from pathlib import Path

import pandas as pd
import pymupdf


ROOT = figure_root()
RAW = private_path("source_pdfs")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fallback-only", action="store_true", help="Refresh only figures marked as raster fallbacks in the manifest")
    args = parser.parse_args()

    fallback_stems = set()
    if args.fallback_only:
        manifest = pd.read_csv(ROOT / "figure_source_manifest.csv")
        fallback_stems = set(manifest.loc[~manifest["status"].astype(str).str.startswith("native vector"), "figure"])

    sources = sorted(RAW.rglob("*.pdf"))
    if not sources:
        raise FileNotFoundError(f"No source PDFs found under {RAW}")

    for source in sources:
        relative = source.relative_to(RAW)
        figure_stem = relative.with_suffix("").as_posix()
        if args.fallback_only and figure_stem not in fallback_stems:
            continue
        output_pdf = ROOT / relative
        output_pdf.parent.mkdir(parents=True, exist_ok=True)
        output_svg = output_pdf.with_suffix("").with_name(output_pdf.stem + ".svg")
        if output_svg.exists() and not args.fallback_only:
            print(f"Keeping existing SVG {output_svg.relative_to(ROOT)}")
            continue
        shutil.copy2(source, output_pdf)

        document = pymupdf.open(source)
        try:
            svg_parts = [page.get_svg_image(text_as_path=False) for page in document]
        finally:
            document.close()

        # Source pages are rasterized; the SVG preserves embedded images rather than editable plot elements.
        output_svg.write_text("\n".join(svg_parts), encoding="utf-8")
        print(f"Wrote {output_pdf.relative_to(ROOT)} and {output_svg.relative_to(ROOT)}")


if __name__ == "__main__":
    main()