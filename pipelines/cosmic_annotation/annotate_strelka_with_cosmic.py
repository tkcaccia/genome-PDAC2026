#!/usr/bin/env python3
import argparse
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def normalize_chrom(chrom: str) -> str:
    return chrom.removeprefix("chr").upper()


def split_info(info: str):
    result = {}
    for token in info.split(";"):
        if "=" in token:
            k, v = token.split("=", 1)
            result[k] = v
        else:
            result[token] = True
    return result


def parse_csq_header(vcf_path: Path):
    with open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ"):
                marker = "Format: "
                idx = line.find(marker)
                if idx == -1:
                    continue
                return line[idx + len(marker):].split('">')[0].split("|")
            if line.startswith("#CHROM"):
                break
    raise RuntimeError(f"CSQ header not found in {vcf_path}")


def extract_first_csq(info_dict, csq_fields):
    raw = info_dict.get("CSQ")
    if not raw:
        return []
    annotations = []
    for entry in raw.split(","):
        vals = entry.split("|")
        if len(vals) < len(csq_fields):
            vals = vals + [""] * (len(csq_fields) - len(vals))
        annotations.append(dict(zip(csq_fields, vals)))
    return annotations


def choose_annotation(annotations):
    if not annotations:
        return {}
    canonical = [a for a in annotations if a.get("CANONICAL") == "YES"]
    if canonical:
        protein = [a for a in canonical if a.get("BIOTYPE") == "protein_coding"]
        return protein[0] if protein else canonical[0]
    protein = [a for a in annotations if a.get("BIOTYPE") == "protein_coding"]
    if protein:
        return protein[0]
    return annotations[0]


def load_cgc(path: Path):
    by_gene = {}
    with open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by_gene[row["GENE_SYMBOL"]] = row
    return by_gene


def load_hallmarks(path: Path):
    by_gene = defaultdict(lambda: {"hallmarks": set(), "descriptions": set()})
    with open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["GENE_SYMBOL"]
            if row["HALLMARK"]:
                by_gene[gene]["hallmarks"].add(row["HALLMARK"])
            if row["DESCRIPTION"]:
                by_gene[gene]["descriptions"].add(row["DESCRIPTION"])
    return by_gene


def load_mutant_census(path: Path):
    by_exact = defaultdict(list)
    with open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (
                row["GENE_SYMBOL"],
                row["HGVSC"],
                row["HGVSP"],
                normalize_chrom(row["CHROMOSOME"]),
                row["GENOME_START"],
                row["GENOMIC_WT_ALLELE"],
                row["GENOMIC_MUT_ALLELE"],
            )
            by_exact[key].append(row)
    return by_exact


def load_resistance(path: Path):
    by_exact = defaultdict(list)
    with open_text(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (
                row["GENE_SYMBOL"],
                row["HGVSC"],
                row["HGVSP"],
                normalize_chrom(row["CHROMOSOME"]),
                row["GENOME_START"],
                row["GENOMIC_WT_ALLELE"],
                row["GENOMIC_MUT_ALLELE"],
            )
            by_exact[key].append(row)
    return by_exact


def summarize_mutant_hits(rows):
    if not rows:
        return {
            "cosmic_mutant_hit": "no",
            "cosmic_mutation_ids": "",
            "cosmic_sample_count": "0",
            "cosmic_mutation_description": "",
            "cosmic_pubmed_pmids": "",
        }
    return {
        "cosmic_mutant_hit": "yes",
        "cosmic_mutation_ids": ";".join(sorted({r["MUTATION_ID"] for r in rows if r["MUTATION_ID"]})),
        "cosmic_sample_count": str(len({r["COSMIC_SAMPLE_ID"] for r in rows if r["COSMIC_SAMPLE_ID"]})),
        "cosmic_mutation_description": ";".join(sorted({r["MUTATION_DESCRIPTION"] for r in rows if r["MUTATION_DESCRIPTION"]})),
        "cosmic_pubmed_pmids": ";".join(sorted({r["PUBMED_PMID"] for r in rows if r["PUBMED_PMID"]})),
    }


def summarize_resistance_hits(rows):
    if not rows:
        return {
            "cosmic_resistance_hit": "no",
            "resistance_drugs": "",
            "resistance_responses": "",
            "resistance_mutation_ids": "",
        }
    return {
        "cosmic_resistance_hit": "yes",
        "resistance_drugs": ";".join(sorted({r["DRUG_NAME"] for r in rows if r["DRUG_NAME"]})),
        "resistance_responses": ";".join(sorted({r["DRUG_RESPONSE"] for r in rows if r["DRUG_RESPONSE"]})),
        "resistance_mutation_ids": ";".join(sorted({r["MUTATION_ID"] for r in rows if r["MUTATION_ID"]})),
    }


def iter_strelka_variants(vcf_path: Path, csq_fields):
    sample_pair = (
        vcf_path.name
        .replace(".strelka.somatic_snvs_VEP.ann.vcf.gz", "")
        .replace(".strelka.somatic_indels_VEP.ann.vcf.gz", "")
    )
    caller_class = "snv" if ".somatic_snvs_" in vcf_path.name else "indel"
    with open_text(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, qual, flt, info, *rest = line.rstrip("\n").split("\t")
            if flt != "PASS":
                continue
            info_dict = split_info(info)
            annotations = extract_first_csq(info_dict, csq_fields)
            chosen = choose_annotation(annotations)
            yield {
                "sample_pair": sample_pair,
                "caller_class": caller_class,
                "chrom": normalize_chrom(chrom),
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "annotation": chosen,
            }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf-dir", required=True)
    parser.add_argument("--cosmic-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    vcf_dir = Path(args.vcf_dir)
    cosmic_dir = Path(args.cosmic_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cgc = load_cgc(cosmic_dir / "Cosmic_CancerGeneCensus_v103_GRCh38.tsv.gz")
    hallmarks = load_hallmarks(cosmic_dir / "Cosmic_CancerGeneCensusHallmarksOfCancer_v103_GRCh38.tsv.gz")
    mutant = load_mutant_census(cosmic_dir / "Cosmic_MutantCensus_v103_GRCh38.tsv.gz")
    resistance = load_resistance(cosmic_dir / "Cosmic_ResistanceMutations_v103_GRCh38.tsv.gz")

    vcfs = sorted(vcf_dir.rglob("*.strelka.somatic_*_VEP.ann.vcf.gz"))
    if not vcfs:
        raise SystemExit(f"No Strelka annotated VCFs found in {vcf_dir}")

    csq_fields = parse_csq_header(vcfs[0])
    rows = []
    gene_counts = Counter()
    exact_counts = Counter()

    for vcf in vcfs:
        for record in iter_strelka_variants(vcf, csq_fields):
            ann = record["annotation"]
            gene = ann.get("SYMBOL", "")
            hgvsc = ann.get("HGVSc", "")
            hgvsp = ann.get("HGVSp", "")
            exact_key = (gene, hgvsc, hgvsp, record["chrom"], record["pos"], record["ref"], record["alt"])

            cgc_row = cgc.get(gene, {})
            hallmark_row = hallmarks.get(gene, {"hallmarks": set(), "descriptions": set()})
            mutant_hit = summarize_mutant_hits(mutant.get(exact_key, []))
            resistance_hit = summarize_resistance_hits(resistance.get(exact_key, []))

            row = {
                "sample_pair": record["sample_pair"],
                "caller_class": record["caller_class"],
                "chrom": record["chrom"],
                "pos": record["pos"],
                "ref": record["ref"],
                "alt": record["alt"],
                "gene_symbol": gene,
                "consequence": ann.get("Consequence", ""),
                "impact": ann.get("IMPACT", ""),
                "canonical": ann.get("CANONICAL", ""),
                "mane_select": ann.get("MANE_SELECT", ""),
                "biotype": ann.get("BIOTYPE", ""),
                "hgvsc": hgvsc,
                "hgvsp": hgvsp,
                "variant_class": ann.get("VARIANT_CLASS", ""),
                "existing_variation": ann.get("Existing_variation", ""),
                "cgc_hit": "yes" if cgc_row else "no",
                "cgc_tier": cgc_row.get("TIER", ""),
                "cgc_role_in_cancer": cgc_row.get("ROLE_IN_CANCER", ""),
                "hallmarks": ";".join(sorted(hallmark_row["hallmarks"])),
                **mutant_hit,
                **resistance_hit,
            }
            rows.append(row)

            if gene:
                gene_counts[gene] += 1
            if mutant_hit["cosmic_mutant_hit"] == "yes":
                exact_counts[(gene, hgvsc, hgvsp)] += 1

    out_main = output_dir / "strelka_pass_cosmic_annotation.tsv"
    with open(out_main, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    out_gene = output_dir / "strelka_pass_recurrent_gene_summary.tsv"
    with open(out_gene, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_symbol", "pass_variant_count"])
        for gene, count in gene_counts.most_common():
            writer.writerow([gene, count])

    out_exact = output_dir / "strelka_pass_exact_cosmic_match_summary.tsv"
    with open(out_exact, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_symbol", "hgvsc", "hgvsp", "match_count"])
        for (gene, hgvsc, hgvsp), count in exact_counts.most_common():
            writer.writerow([gene, hgvsc, hgvsp, count])


if __name__ == "__main__":
    main()
