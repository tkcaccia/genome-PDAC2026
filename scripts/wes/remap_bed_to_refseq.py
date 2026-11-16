#!/usr/bin/env python3
# WES processing, variant summaries, coverage and paired MSI.
# Private paths/cohorts: config/private.example.json; never store populated inputs here.
"""Publication analysis; inputs and outputs are documented in docs/workflow.md."""
from pdac_runtime import project_root, private_path, private_value, cohort as configured_cohort, code_path
import re
import sys


def build_contig_map(dict_path):
    contig_map = {}
    with open(dict_path) as handle:
        for line in handle:
            if not line.startswith("@SQ"):
                continue
            name = next((field[3:] for field in line.rstrip().split("\t") if field.startswith("SN:")), None)
            if name is None:
                continue
            match = re.fullmatch(r"NC_0*(\d+)\.\d+", name)
            if not match:
                continue
            number = int(match.group(1))
            if 1 <= number <= 22:
                contig_map[f"chr{number}"] = name
                contig_map[str(number)] = name
            elif number == 23:
                contig_map.update({"chrX": name, "X": name})
            elif number == 24:
                contig_map.update({"chrY": name, "Y": name})
            elif number == 12920:
                contig_map.update({"chrM": name, "MT": name, "M": name})
    return contig_map


def remap_bed(dict_path, bed_in, bed_out):
    contig_map = build_contig_map(dict_path)
    if not contig_map:
        raise ValueError(f"No primary GRCh38 RefSeq contigs found in {dict_path}")

    kept = 0
    with open(bed_in) as source, open(bed_out, "w") as target:
        for line_number, line in enumerate(source, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed BED row {line_number}: expected at least 3 columns")
            chromosome = fields[0]
            mapped = contig_map.get(chromosome)
            if mapped is None and chromosome.startswith("chr"):
                mapped = contig_map.get(chromosome[3:])
            if mapped is None:
                raise ValueError(f"Unmapped contig {chromosome!r} at BED row {line_number}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid coordinates at BED row {line_number}: {start}-{end}")
            fields[0] = mapped
            target.write("\t".join(fields) + "\n")
            kept += 1

    if kept == 0:
        raise ValueError("No BED intervals were written")
    print(f"Remapped {kept} intervals to RefSeq contig names: {bed_out}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("Usage: remap_bed_to_refseq.py reference.dict input.bed output.bed")
    remap_bed(*sys.argv[1:])