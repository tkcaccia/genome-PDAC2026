#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge, PathPatch, Patch
from matplotlib.path import Path as MplPath


ALT_BREAKEND_RE = re.compile(r"[\[\]]([^:\[\]]+):([0-9]+)[\[\]]")
CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PDAC_GENES = [
    "KRAS",
    "TP53",
    "CDKN2A",
    "SMAD4",
    "ARID1A",
    "RNF43",
    "BRCA1",
    "BRCA2",
    "PALB2",
    "ATM",
    "STK11",
    "GNAS",
    "TGFBR2",
    "KDM6A",
    "KMT2C",
    "MAP2K4",
    "PBRM1",
    "BRAF",
    "ERBB2",
    "MYC",
    "CCND1",
]
HIGH_IMPACT_CONSEQUENCES = {
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate per-patient circos-style plots combining SNV, CNV, and SV calls."
    )
    parser.add_argument("--mutect2-root", required=True)
    parser.add_argument("--ascat-root", required=True)
    parser.add_argument("--manta-root", required=True)
    parser.add_argument("--vep-root")
    parser.add_argument("--star-root")
    parser.add_argument("--gtf")
    parser.add_argument("--cytoband")
    parser.add_argument("--fai", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--bin-size", type=int, default=5_000_000)
    parser.add_argument("--samples", nargs="*", default=[])
    return parser.parse_args()


def read_fai(path: Path):
    chrom_lengths = {}
    with path.open() as fh:
        for line in fh:
            chrom, length, *_ = line.rstrip().split("\t")
            if chrom in CHROM_ORDER:
                chrom_lengths[chrom] = int(length)
    return {chrom: chrom_lengths[chrom] for chrom in CHROM_ORDER if chrom in chrom_lengths}


def build_genome_layout(chrom_lengths, gap_deg=2.0):
    total_len = sum(chrom_lengths.values())
    usable_deg = 360.0 - gap_deg * len(chrom_lengths)
    layout = {}
    cur_deg = 90.0
    for chrom, length in chrom_lengths.items():
        span = usable_deg * (length / total_len)
        start_deg = cur_deg
        end_deg = cur_deg + span
        layout[chrom] = {
            "length": length,
            "start": math.radians(start_deg),
            "end": math.radians(end_deg),
            "start_deg": start_deg,
            "end_deg": end_deg,
            "span_deg": span,
        }
        cur_deg += span + gap_deg
    return layout


def chrom_pos_to_theta(layout, chrom, pos):
    info = layout[chrom]
    frac = min(max(pos, 1), info["length"]) / info["length"]
    return info["start"] + (info["end"] - info["start"]) * frac


def polar_to_xy(theta, radius):
    return radius * math.cos(theta), radius * math.sin(theta)


def parse_info(info_str):
    info = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            info[k] = v
        else:
            info[field] = True
    return info


def parse_gtf_attrs(attr_str):
    attrs = {}
    for field in attr_str.rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue
        if " " not in field:
            continue
        key, value = field.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def read_gene_coords(gtf_path, chrom_lengths, gene_names):
    if not gtf_path:
        return {}
    coords = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, _source, feature, start, end, _score, _strand, _phase, attrs = line.rstrip().split("\t")
            if feature != "gene" or chrom not in chrom_lengths:
                continue
            parsed = parse_gtf_attrs(attrs)
            gene_name = parsed.get("gene_name")
            if gene_name not in gene_names:
                continue
            coords[gene_name] = {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "gene_id": parsed.get("gene_id", "").split(".")[0],
            }
    return coords


def read_cytobands(path, chrom_lengths):
    if not path:
        return defaultdict(list)
    bands = defaultdict(list)
    with Path(path).open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, name, stain = line.rstrip().split("\t")[:5]
            if chrom not in chrom_lengths:
                continue
            bands[chrom].append((int(start) + 1, int(end), name, stain))
    return bands


def cytoband_color(stain):
    return {
        "gneg": "#ffffff",
        "gpos25": "#d9d9d9",
        "gpos50": "#a6a6a6",
        "gpos75": "#737373",
        "gpos100": "#252525",
        "gvar": "#e6e6e6",
        "stalk": "#bdbdbd",
        "acen": "#e34a33",
    }.get(stain, "#f0f0f0")


def patient_id_from_pair(sample_pair):
    match = re.match(r"(.+?)T_WES_vs_.+?N_WES$", sample_pair)
    return match.group(1) if match else sample_pair.split("T_WES", 1)[0]


def read_star_counts(star_root, sample_id):
    if not star_root:
        return {}, 0.0
    path = Path(star_root) / f"{sample_id}.ReadsPerGene.out.tab"
    if not path.exists():
        return {}, 0.0
    counts = {}
    library_total = 0.0
    with path.open() as fh:
        for line in fh:
            gene_id, unstranded, *_rest = line.rstrip().split("\t")
            if gene_id.startswith("N_"):
                continue
            value = float(unstranded)
            counts[gene_id.split(".")[0]] = value
            library_total += value
    return counts, library_total


def expression_log2fc_for_pair(star_root, sample_pair, gene_coords):
    if not star_root:
        return {}
    pid = patient_id_from_pair(sample_pair)
    tumor_counts, tumor_total = read_star_counts(star_root, f"{pid}T_RNA")
    normal_counts, normal_total = read_star_counts(star_root, f"{pid}N_RNA")
    if tumor_total == 0 or normal_total == 0:
        return {}
    result = {}
    for gene, meta in gene_coords.items():
        gene_id = meta.get("gene_id")
        if not gene_id:
            continue
        tumor_cpm = (tumor_counts.get(gene_id, 0.0) + 0.5) / tumor_total * 1_000_000
        normal_cpm = (normal_counts.get(gene_id, 0.0) + 0.5) / normal_total * 1_000_000
        result[gene] = math.log2(tumor_cpm / normal_cpm)
    return result


def expr_color(log2fc):
    if log2fc is None:
        return "#969696"
    if log2fc >= 1:
        return "#b2182b"
    if log2fc <= -1:
        return "#2166ac"
    return "#f7f7f7"


def read_vep_gene_hits(vcf_path, genes):
    hits = defaultdict(list)
    if not vcf_path or not Path(vcf_path).exists():
        return hits
    csq_fields = []
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                fmt = line.split("Format: ", 1)[1].split('">', 1)[0]
                csq_fields = fmt.split("|")
                continue
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, qual, filt, info_str, *_rest = line.rstrip().split("\t")
            if filt != "PASS":
                continue
            if not csq_fields:
                continue
            info = parse_info(info_str)
            for annot in str(info.get("CSQ", "")).split(","):
                values = annot.split("|")
                csq = {field: values[i] if i < len(values) else "" for i, field in enumerate(csq_fields)}
                symbol = csq.get("SYMBOL", "")
                if symbol not in genes:
                    continue
                consequences = set(csq.get("Consequence", "").split("&"))
                impact = csq.get("IMPACT", "")
                is_disruptive = impact == "HIGH" or bool(consequences & HIGH_IMPACT_CONSEQUENCES)
                hits[symbol].append(
                    {
                        "chrom": chrom,
                        "pos": int(pos),
                        "impact": impact,
                        "consequence": csq.get("Consequence", ""),
                        "disruptive": is_disruptive,
                    }
                )
    return hits


def parse_alt_breakend(alt):
    m = ALT_BREAKEND_RE.search(alt)
    if not m:
        return None, None
    return m.group(1), int(m.group(2))


def parse_mutect2_snvs(vcf_path, chrom_lengths, bin_size):
    density = defaultdict(Counter)
    total = 0
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, qual, filt, *_rest = line.rstrip().split("\t")
            if chrom not in chrom_lengths:
                continue
            if filt != "PASS":
                continue
            alts = alt.split(",")
            if len(ref) != 1:
                continue
            if any(len(a) != 1 or a not in "ACGT" for a in alts):
                continue
            pos = int(pos)
            density[chrom][(pos - 1) // bin_size] += 1
            total += 1
    return density, total


def parse_ascat_segments(seg_path, chrom_lengths):
    rows = []
    with seg_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chrom = row["chr"]
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            if chrom not in chrom_lengths:
                continue
            start = int(float(row["startpos"]))
            end = int(float(row["endpos"]))
            nmajor = int(round(float(row["nMajor"])))
            nminor = int(round(float(row["nMinor"])))
            rows.append((chrom, start, end, nmajor, nminor))
    return rows


def parse_manta_events(vcf_path, chrom_lengths):
    events = []
    seen = set()
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom, pos, record_id, ref, alt, qual, filt, info_str, *_rest = line.rstrip().split("\t")
            if chrom not in chrom_lengths:
                continue
            info = parse_info(info_str)
            svtype = info.get("SVTYPE", "")
            pos = int(pos)

            if svtype == "BND":
                mate_chrom, mate_pos = parse_alt_breakend(alt)
                if mate_chrom not in chrom_lengths or mate_pos is None:
                    continue
                mate_id = info.get("MATEID", "")
                key = tuple(sorted([record_id, mate_id])) if mate_id else (record_id,)
                if key in seen:
                    continue
                seen.add(key)
                events.append(
                    {
                        "chrom1": chrom,
                        "pos1": pos,
                        "chrom2": mate_chrom,
                        "pos2": mate_pos,
                        "svtype": "TRA" if mate_chrom != chrom else "BND",
                        "filter": filt,
                    }
                )
                continue

            if "END" not in info:
                continue
            end = int(info["END"])
            if chrom not in chrom_lengths:
                continue
            events.append(
                {
                    "chrom1": chrom,
                    "pos1": pos,
                    "chrom2": chrom,
                    "pos2": end,
                    "svtype": svtype,
                    "filter": filt,
                }
            )
    return events


def cn_color(total_cn):
    if total_cn <= 0:
        return "#08306b"
    if total_cn == 1:
        return "#6baed6"
    if total_cn == 2:
        return "#d9d9d9"
    if total_cn == 3:
        return "#fdae6b"
    return "#cb181d"


def segment_overlaps_gene(segment, gene):
    chrom, start, end, nmajor, nminor = segment
    return chrom == gene["chrom"] and start <= gene["end"] and end >= gene["start"]


def add_arc(ax, theta1, theta2, radius_outer, radius_inner, color, lw=0.6, alpha=0.8):
    x1, y1 = polar_to_xy(theta1, radius_outer)
    x2, y2 = polar_to_xy(theta2, radius_outer)
    c1x, c1y = polar_to_xy(theta1, radius_inner)
    c2x, c2y = polar_to_xy(theta2, radius_inner)
    verts = [(x1, y1), (c1x, c1y), (c2x, c2y), (x2, y2)]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    patch = PathPatch(MplPath(verts, codes), facecolor="none", edgecolor=color, lw=lw, alpha=alpha)
    ax.add_patch(patch)


def draw_sample(
    sample_pair,
    mutect2_vcf,
    ascat_seg,
    manta_vcf,
    vep_vcf,
    chrom_lengths,
    layout,
    output_path,
    bin_size,
    cytobands,
    gene_coords,
    star_root,
):
    snv_density, snv_total = parse_mutect2_snvs(mutect2_vcf, chrom_lengths, bin_size)
    snv_by_chrom = {chrom: sum(snv_density.get(chrom, {}).values()) for chrom in chrom_lengths}
    segments = parse_ascat_segments(ascat_seg, chrom_lengths)
    sv_events = parse_manta_events(manta_vcf, chrom_lengths)
    gene_hits = read_vep_gene_hits(vep_vcf, set(gene_coords))
    gene_log2fc = expression_log2fc_for_pair(star_root, sample_pair, gene_coords)
    translocations = [e for e in sv_events if e["chrom1"] != e["chrom2"]]
    intrachrom_sv = [e for e in sv_events if e["chrom1"] == e["chrom2"]]
    biallelic_genes = set()
    for gene, coord in gene_coords.items():
        has_disruptive_variant = any(hit["disruptive"] for hit in gene_hits.get(gene, []))
        has_loss = any(segment_overlaps_gene(seg, coord) and (seg[4] == 0 or seg[3] + seg[4] <= 1) for seg in segments)
        if has_disruptive_variant and has_loss:
            biallelic_genes.add(gene)

    fig, ax = plt.subplots(figsize=(13, 13))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-1.42, 1.42)
    ax.set_ylim(-1.42, 1.42)

    r_outer = 1.02
    ideogram_width = 0.08
    gene_r0, gene_r1 = 0.93, 1.00
    load_r0, load_r1 = 0.82, 0.92
    snv_r0, snv_r1 = 0.68, 0.80
    cnv_r0, cnv_r1 = 0.50, 0.66
    link_r = 0.46
    label_r = 1.15

    for i, (chrom, info) in enumerate(layout.items()):
        if cytobands.get(chrom):
            for start, end, band_name, stain in cytobands[chrom]:
                th1 = math.degrees(chrom_pos_to_theta(layout, chrom, start))
                th2 = math.degrees(chrom_pos_to_theta(layout, chrom, end))
                ax.add_patch(
                    Wedge(
                        (0, 0),
                        r_outer,
                        th1,
                        th2,
                        width=ideogram_width,
                        facecolor=cytoband_color(stain),
                        edgecolor="#f2f2f2",
                        lw=0.15,
                    )
                )
        else:
            face = "#f0f0f0" if i % 2 == 0 else "#d9d9d9"
            ax.add_patch(
                Wedge(
                    (0, 0),
                    r_outer,
                    info["start_deg"],
                    info["end_deg"],
                    width=ideogram_width,
                    facecolor=face,
                    edgecolor="white",
                    lw=1,
                )
            )
        mid = (info["start"] + info["end"]) / 2
        lx, ly = polar_to_xy(mid, label_r)
        ax.text(lx, ly, chrom.replace("chr", ""), ha="center", va="center", fontsize=8)

        # Track backgrounds make chromosomes with zero calls visible too.
        ax.add_patch(
            Wedge(
                (0, 0),
                load_r1,
                info["start_deg"],
                info["end_deg"],
                width=(load_r1 - load_r0),
                facecolor="#f7f7f7",
                edgecolor="white",
                lw=0.25,
            )
        )
        ax.add_patch(
            Wedge(
                (0, 0),
                gene_r1,
                info["start_deg"],
                info["end_deg"],
                width=(gene_r1 - gene_r0),
                facecolor="#f7f7f7",
                edgecolor="white",
                lw=0.25,
            )
        )
        ax.add_patch(
            Wedge(
                (0, 0),
                snv_r1,
                info["start_deg"],
                info["end_deg"],
                width=(snv_r1 - snv_r0),
                facecolor="#fff5eb",
                edgecolor="white",
                lw=0.25,
            )
        )

    for gene in PDAC_GENES:
        coord = gene_coords.get(gene)
        if not coord:
            continue
        theta = chrom_pos_to_theta(layout, coord["chrom"], (coord["start"] + coord["end"]) // 2)
        log2fc = gene_log2fc.get(gene)
        x0, y0 = polar_to_xy(theta, gene_r0 + 0.012)
        x1, y1 = polar_to_xy(theta, gene_r1 - 0.012)
        color = expr_color(log2fc)
        ax.plot([x0, x1], [y0, y1], color=color, lw=3.0, solid_capstyle="round", zorder=8)
        if gene in biallelic_genes:
            sx, sy = polar_to_xy(theta, gene_r1 + 0.025)
            ax.scatter([sx], [sy], marker="*", s=55, color="#000000", zorder=9)
        if gene in biallelic_genes or gene in gene_hits or (log2fc is not None and abs(log2fc) >= 1):
            tx, ty = polar_to_xy(theta, 1.25)
            ha = "left" if math.cos(theta) >= 0 else "right"
            ax.text(tx, ty, gene, ha=ha, va="center", fontsize=7, color="#252525")

    max_chrom_load = max(snv_by_chrom.values(), default=1) or 1
    for chrom, count in snv_by_chrom.items():
        if count == 0:
            continue
        info = layout[chrom]
        height = (count / max_chrom_load) * (load_r1 - load_r0)
        ax.add_patch(
            Wedge(
                (0, 0),
                load_r0 + height,
                info["start_deg"],
                info["end_deg"],
                width=height,
                facecolor="#d94801",
                edgecolor="none",
                alpha=0.95,
            )
        )

    max_bin = max((count for bins in snv_density.values() for count in bins.values()), default=1)
    for chrom, bins in snv_density.items():
        chrom_len = chrom_lengths[chrom]
        nbins = math.ceil(chrom_len / bin_size)
        for idx in range(nbins):
            count = bins.get(idx, 0)
            if count == 0:
                continue
            start = idx * bin_size + 1
            end = min((idx + 1) * bin_size, chrom_len)
            th1 = math.degrees(chrom_pos_to_theta(layout, chrom, start))
            th2 = math.degrees(chrom_pos_to_theta(layout, chrom, end))
            height = (count / max_bin) * (snv_r1 - snv_r0)
            ax.add_patch(Wedge((0, 0), snv_r0 + height, th1, th2, width=height, facecolor="#8c2d04", edgecolor="none", alpha=0.9))

    for chrom, start, end, nmajor, nminor in segments:
        total_cn = nmajor + nminor
        th1 = math.degrees(chrom_pos_to_theta(layout, chrom, start))
        th2 = math.degrees(chrom_pos_to_theta(layout, chrom, end))
        ax.add_patch(Wedge((0, 0), cnv_r1, th1, th2, width=(cnv_r1 - cnv_r0), facecolor=cn_color(total_cn), edgecolor="none", alpha=0.95))

    for event in intrachrom_sv:
        theta1 = chrom_pos_to_theta(layout, event["chrom1"], event["pos1"])
        theta2 = chrom_pos_to_theta(layout, event["chrom2"], event["pos2"])
        color = "#3182bd" if event["svtype"] != "BND" else "#756bb1"
        add_arc(ax, theta1, theta2, link_r, 0.18, color, lw=0.7, alpha=0.45)

    for event in translocations:
        theta1 = chrom_pos_to_theta(layout, event["chrom1"], event["pos1"])
        theta2 = chrom_pos_to_theta(layout, event["chrom2"], event["pos2"])
        add_arc(ax, theta1, theta2, link_r, 0.08, "#e41a1c", lw=0.9, alpha=0.8)

    title = (
        f"{sample_pair}\n"
        f"SNVs(PASS): {snv_total} | CNV segments: {len(segments)} | "
        f"SV events: {len(sv_events)} | interchromosomal translocations: {len(translocations)} | "
        f"biallelic PDAC genes: {len(biallelic_genes)}"
    )
    ax.text(0, 1.30, title, ha="center", va="center", fontsize=12, fontweight="bold")
    ax.text(
        0,
        0,
        f"PASS SNV load\n{snv_total:,} SNVs\nmax {max_bin:,}/{bin_size // 1_000_000} Mb bin",
        ha="center",
        va="center",
        fontsize=12,
        fontweight="bold",
        color="#4a1d05",
    )

    legend_handles = [
        Line2D([0], [0], color="#d94801", lw=6, label="Chromosome SNV load"),
        Line2D([0], [0], color="#8c2d04", lw=6, label=f"SNV density ({bin_size // 1_000_000} Mb bins)"),
        Line2D([0], [0], color="#b2182b", lw=3, label="PDAC gene RNA up"),
        Line2D([0], [0], color="#2166ac", lw=3, label="PDAC gene RNA down"),
        Line2D([0], [0], marker="*", color="#000000", lw=0, markersize=8, label="Likely biallelic gene"),
        Line2D([0], [0], color="#3182bd", lw=2, label="Intrachromosomal SV"),
        Line2D([0], [0], color="#e41a1c", lw=2, label="Interchromosomal translocation"),
        Line2D([0], [0], color="#756bb1", lw=2, label="Intrachromosomal BND"),
        Patch(facecolor="#e34a33", edgecolor="none", label="Centromere/cytoband"),
        Line2D([0], [0], color=cn_color(1), lw=6, label="CN loss"),
        Line2D([0], [0], color=cn_color(2), lw=6, label="CN neutral"),
        Line2D([0], [0], color=cn_color(3), lw=6, label="CN gain"),
        Line2D([0], [0], color=cn_color(4), lw=6, label="CN amplification"),
    ]
    ax.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.5, -0.045), ncol=4, frameon=False, fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    chrom_lengths = read_fai(Path(args.fai))
    layout = build_genome_layout(chrom_lengths)

    mutect2_root = Path(args.mutect2_root)
    ascat_root = Path(args.ascat_root)
    manta_root = Path(args.manta_root)
    vep_root = Path(args.vep_root) if args.vep_root else None
    star_root = Path(args.star_root) if args.star_root else None
    cytobands = read_cytobands(args.cytoband, chrom_lengths)
    gene_coords = read_gene_coords(Path(args.gtf), chrom_lengths, set(PDAC_GENES)) if args.gtf else {}

    mutect2_map = {
        p.name.replace(".mutect2.filtered.vcf.gz", ""): p
        for p in mutect2_root.rglob("*.mutect2.filtered.vcf.gz")
    }
    ascat_map = {
        p.name.replace(".segments.txt", ""): p
        for p in ascat_root.rglob("*.segments.txt")
    }
    manta_map = {
        p.name.replace(".manta.somatic_sv.vcf.gz", ""): p
        for p in manta_root.rglob("*.manta.somatic_sv.vcf.gz")
    }
    vep_map = {}
    if vep_root:
        vep_map = {
            p.name.replace(".mutect2.filtered_VEP.ann.vcf.gz", ""): p
            for p in vep_root.rglob("*.mutect2.filtered_VEP.ann.vcf.gz")
        }

    sample_pairs = sorted(set(mutect2_map) & set(ascat_map) & set(manta_map))
    if args.samples:
        requested = set(args.samples)
        sample_pairs = [s for s in sample_pairs if s in requested]

    manifest_rows = []
    for sample_pair in sample_pairs:
        out_png = output_dir / f"{sample_pair}.circos.png"
        draw_sample(
            sample_pair=sample_pair,
            mutect2_vcf=mutect2_map[sample_pair],
            ascat_seg=ascat_map[sample_pair],
            manta_vcf=manta_map[sample_pair],
            vep_vcf=vep_map.get(sample_pair),
            chrom_lengths=chrom_lengths,
            layout=layout,
            output_path=out_png,
            bin_size=args.bin_size,
            cytobands=cytobands,
            gene_coords=gene_coords,
            star_root=star_root,
        )
        manifest_rows.append({"sample_pair": sample_pair, "plot_png": str(out_png)})

    with (output_dir / "manifest.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_pair", "plot_png"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)


if __name__ == "__main__":
    main()
