# Reviewer Comments for Journal of Translational Medicine

## Overall Assessment

This manuscript addresses an important and underrepresented area: integrated molecular profiling of pancreatic cancer in African patients. The strongest aspect is the multi-layer integration of WES, RNA-seq, copy-number, structural-variant, fusion, immune/stromal, KRAS and MSI/MMR analyses into a patient-level interpretation framework. The potential identification of a KRAS-wild-type MSI/MMR-deficient candidate is biologically and translationally interesting.

However, the manuscript remains exploratory. The cohort is small, several findings are computational only, and the claims must be calibrated carefully. The paper could be suitable for Journal of Translational Medicine if it is framed as a feasibility/hypothesis-generating multi-omics study and if key methods and validation limitations are made fully transparent.

## Major Comments

1. The KRAS-wild-type claim needs very careful wording. The lower KRAS alteration frequency is interesting, but in a 14-patient cohort it should not be presented as an African-population difference. The authors should emphasize sample size, tumour purity, histological confirmation, coverage, FFPE/calling artifacts and the need for orthogonal validation.

2. Patient 23 is the strongest translational finding. The manuscript should elevate this as a single high-priority computational candidate rather than making broad cohort-level claims. MMR IHC for MLH1, PMS2, MSH2 and MSH6, plus MSI-PCR or an orthogonal clinical MSI assay, would be the ideal validation if tissue becomes available.

3. The hypermutation narrative should be softened. The current draft appropriately notes that broad hypermutation was attenuated by stricter filtering. This should remain a limitation, not a central conclusion.

4. Differential-expression methods require explicit clarification. If any reported table came from a custom Python log2CPM/t-test fallback, it should not be called DESeq2, edgeR or limma. The manuscript should state exactly which final expression results were generated with limma-voom or should move the fallback outputs to exploratory/supplementary status.

5. Immune deconvolution should be interpreted conservatively. ESTIMATE, MCP-counter, CIBERSORT, EPIC, xCell and quanTIseq use different scales and assumptions. The manuscript should avoid treating all outputs as percentages and should separate nominal trends from FDR-significant findings.

6. RNA fusion results need clearer completeness language. The manuscript should state which nf-core/rnafusion modules completed, which failed or were skipped, and how candidate fusions were filtered.

7. Copy-number and SV conclusions should distinguish WES-compatible inference from WGS-grade inference. WES-derived CNV/SV analyses are useful but have lower resolution for genome-wide structural patterns than whole-genome sequencing.

8. The actionable-alteration section should avoid treatment recommendations. Phrases such as "potentially relevant", "exploratory" and "may warrant validation" are appropriate.

9. The African cohort framing needs context but not overreach. The underrepresentation of African patients is important, but the study does not have a matched European cohort processed identically. Any comparison to European-heavy cohorts should be presented as literature-contextual, not a formal ancestry association.

10. The methods should include a transparent data-flow diagram or table. Readers need to know which files came from nf-core/rnaseq, nf-core/sarek, nf-core/rnafusion, SigProfilerAssignment, MSIsensor-pro and custom scripts.

## Minor Comments

1. Add exact versions for all downstream R packages where possible.

2. Define abbreviations at first use: TME, CAF, EMT, MSI, MMR, TMB, CNA, SV and DDR.

3. State whether samples are all PDAC histology or broader pancreatic cancer cases.

4. Clarify whether "tumour" and "normal" labels derive from sample names, metadata or pathology reports.

5. Add a concise table separating supported, candidate and unvalidated findings.

6. Make the figure legends self-contained, especially for deconvolution scores that are not percentages.

7. Include the exact TMB denominator used for WES estimates.

8. Clarify whether germline variants are clinically confirmed or research-grade annotations.

9. Include a statement that GitHub contains code and methods only, not patient data.

10. If manuscript figures include patient identifiers, consider anonymized patient labels consistently.

## Recommendation

My reviewer recommendation would be major revision before acceptance. The study is promising and clinically interesting, especially around patient 23, but the manuscript should be tightened around a few defensible claims: noncanonical KRAS-wild-type candidates, one computational MSI/MMR-deficient case, and heterogeneous immune/stromal phenotypes in an underrepresented cohort.
