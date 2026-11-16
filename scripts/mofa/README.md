# MOFA Analysis

One model row represents one primary-cohort patient. The three views are a
retained-gene binary WES matrix, paired RNA log2 responses and prepared tumour-only
protein log2 abundances. Matched-normal protein abundance is not an additional
MOFA specimen or a protein response in this fit.

| Source | Role in the manuscript |
| --- | --- |
| `run_mofa_revised_11.py` | Align views, prepare QC matrices, fit MOFA, export scores/loadings and select clusters; Figure 5B-C |
| `formal_de_mutation_overlap.py` | Intersect formal paired-DESeq2 genes with mutated genes; Figure 5A, not an expression-threshold substitute |
| `pathways_and_targets.py` | Factor-associated RNA enrichment and identifier-mapped cross-modal ranking |
| `finalize_candidates.py` | Candidate/pathway selection supporting factor and candidate interpretation |
| `regenerate_figure11_composite_score.py` | Percentile-normalised, modality-supported candidate scores; Figure 5D |
| `manuscript_clustering.py` | Cluster scatter presentation; Figure 5C |
| `generate_mofa_seed_stability_11.py` | Numerical stability against the saved primary model on identical feature/sample names and QC values |
| `model_inputs.py` | Read-only saved-model layout, axis validation and matrix alignment; no model fitting or feature reselection |
| `generate_figure_s10_reconstruction.py` | Figure S10: stored-name alignment, per-factor and full-model training reconstruction |
| `generate_remaining_mofa_publication_figures_11.py` | Reconstruction and top-200 modality-support summary; Figure S12 |
| `export_modality_support_input.py` | Exact top-200 numeric input for Figure S12, without the unreported top-200 scatter/pie panel |
| `cluster_stability_panel.py` | Fixed-factor sample-exclusion heatmap and silhouette summary; Figure S11 |
| `generate_kras_integrated_figures_11.py` | Integrated feature heatmap and KRAS-stratified feature tables; Figure S13 and corresponding results |

Fitting, score export and plotting are different operations. Reconstruction uses
the saved model and its exact sample/feature ordering; it does not refit MOFA.
Exclusion analysis reclusters retained rows of existing factor scores, not a
retrained model. Stability inputs are taken from the saved model's names and the
corresponding prepared QC matrices; a second feature-selection procedure is not
performed. Correlations are calculated against saved primary-model factors, not
only against another sensitivity fit.

The stability implementation writes to `mofa/seed_stability_model_aligned/` and
refuses an existing directory. It does not overwrite the supplied stability
tables, model or publication figures. Its use requires an explicit scientific
rerun before any manuscript stability values are replaced. Synthetic alignment
tests establish input-handling behaviour, not biological or numerical validation
of the study's fitted model.

The public folder contains code only. Keep the model, molecular matrices, source
tables, SVG/PDF figures, exclusions and real sample ordering outside the repository.
Source hashes and configuration adaptations are recorded in the
[code manifest](../../docs/source_manifest.json). Extra top-100/top-200 plotting
scripts, response-threshold Venn plots, t-SNE and report-writing utilities are not
included as standalone publication analyses.
