#!/usr/bin/env python3
"""Report basic manuscript figure dimensions without reading patient tables.

The script inspects raster image headers and writes a TSV with pixel dimensions,
approximate print size at a target DPI, and simple readiness flags. It is meant
for reproducibility/documentation; it does not upload or modify figures.
"""

from __future__ import annotations

import argparse
import struct
from pathlib import Path


def png_size(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        raise ValueError(f"{path} is not a PNG file")
    width, height = struct.unpack(">II", header[16:24])
    return int(width), int(height)


def image_size(path: Path) -> tuple[int, int]:
    if path.suffix.lower() == ".png":
        return png_size(path)
    raise ValueError(f"Unsupported image type for {path}: {path.suffix}")


def qc_status(width: int, height: int, min_short_edge: int) -> str:
    return "pass" if min(width, height) >= min_short_edge else "review"


def main() -> int:
    parser = argparse.ArgumentParser(description="Create a simple figure QC TSV.")
    parser.add_argument("figure_dir", type=Path, help="Directory containing PNG figures")
    parser.add_argument("output_tsv", type=Path, help="Output TSV path")
    parser.add_argument("--target-dpi", type=float, default=300.0, help="Target print DPI")
    parser.add_argument(
        "--min-short-edge",
        type=int,
        default=1200,
        help="Minimum short-edge pixels for a simple pass/review flag",
    )
    args = parser.parse_args()

    figures = sorted(args.figure_dir.glob("*.png"))
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with args.output_tsv.open("w", encoding="utf-8") as out:
        out.write(
            "\t".join(
                [
                    "figure",
                    "width_px",
                    "height_px",
                    "short_edge_px",
                    "long_edge_px",
                    "print_width_in_at_target_dpi",
                    "print_height_in_at_target_dpi",
                    "target_dpi",
                    "qc_status",
                    "notes",
                ]
            )
            + "\n"
        )
        for figure in figures:
            width, height = image_size(figure)
            short_edge = min(width, height)
            long_edge = max(width, height)
            status = qc_status(width, height, args.min_short_edge)
            notes = (
                "short edge meets simple raster-size screen"
                if status == "pass"
                else "review resolution or export at larger size"
            )
            out.write(
                "\t".join(
                    [
                        figure.name,
                        str(width),
                        str(height),
                        str(short_edge),
                        str(long_edge),
                        f"{width / args.target_dpi:.2f}",
                        f"{height / args.target_dpi:.2f}",
                        f"{args.target_dpi:.0f}",
                        status,
                        notes,
                    ]
                )
                + "\n"
            )

    print(f"Wrote {args.output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
