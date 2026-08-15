# Figure QC

Reusable manuscript-figure quality-control helper.

`figure_qc.py` inspects PNG figure headers and writes a TSV with:

- pixel dimensions
- short and long edge size
- approximate print size at a target DPI, default `300`
- simple `pass`/`review` flag based on minimum short-edge pixels

Example:

```bash
python3 pipelines/figure_qc/figure_qc.py \
  path/to/figures \
  path/to/figure_qc.tsv
```

The script reports dimensions only. It does not upload, modify or embed patient-derived figures.
