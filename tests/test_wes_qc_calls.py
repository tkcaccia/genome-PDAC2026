import gzip
import importlib.util
from pathlib import Path


MODULE_PATH = (
    Path(__file__).parents[1]
    / "pipelines"
    / "wes_qc_calls"
    / "summarize_wes_qc_and_calls.py"
)
SPEC = importlib.util.spec_from_file_location("wes_qc_calls", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
import sys
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_interval_merge_intersection_and_lookup():
    left = {"chr1": MODULE.merge_intervals([(10, 20), (18, 25), (40, 50)])}
    right = {"chr1": MODULE.merge_intervals([(15, 22), (45, 60)])}
    overlap = MODULE.intersect_interval_sets(left, right)
    assert overlap == {"chr1": [(15, 22), (45, 50)]}
    assert MODULE.interval_bases(overlap) == 12
    index = MODULE.IntervalIndex(overlap)
    assert index.contains("chr1", 15)
    assert index.contains("chr1", 49)
    assert not index.contains("chr1", 22)


def test_vcf_summary_restricts_to_pass_and_callable(tmp_path):
    vcf = tmp_path / "test.vcf.gz"
    content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENT-N\tPATIENT-T
chr1\t11\t.\tA\tG\t100\tPASS\tIMPACT=TP53,TX,missense_variant,MISSENSE;TIER=HOTSPOT\tGT:AD:AF:DP\t0/0:30,0:0:30\t0/1:20,10:0.333:30
chr1\t21\t.\tA\tAT\t50\tminTumorVAF\tIMPACT=GENE,TX,frameshift_variant,NONSENSE_OR_FRAMESHIFT\tGT:AD:AF:DP\t0/0:30,0:0:30\t0/1:29,1:0.033:30
chr2\t11\t.\tC\tT\t80\tPASS\tIMPACT=GENE2,TX,intron_variant,NONE\tGT:AD:AF:DP\t0/0:30,0:0:30\t0/1:20,10:0.333:30
"""
    with gzip.open(vcf, "wt") as handle:
        handle.write(content)
    index = MODULE.IntervalIndex({"chr1": [(0, 20)]})
    summary, filters, calls = MODULE.summarize_vcf(vcf, "1", "PURPLE", index, True)
    assert summary["vcf_records"] == 3
    assert summary["pass_pair_callable_records"] == 1
    assert summary["pass_pair_callable_nonsynonymous"] == 1
    assert summary["pass_pair_callable_hotspots"] == 1
    assert summary["strict_read_evidence_records"] == 1
    assert summary["strict_read_evidence_nonsynonymous"] == 1
    assert filters["PASS"] == 2
    assert calls[0]["gene"] == "TP53"
    assert calls[0]["tumour_alt_depth"] == 10
    assert calls[0]["strict_read_evidence"] is True
