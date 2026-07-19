# Analysis Code — Supplementary Materials

## Study Title
Systematic Mendelian Randomization Analysis of Environmental Exposures and Chronic Diseases: Evidence from Univariate MR, Multivariable MR, Mediation MR, and Colocalization

## Requirements
- R >= 4.5.0
- Python >= 3.10
- PLINK 1.9 (for LD clumping)
- R packages: TwoSampleMR, ieugwasr, coloc, MRPRESSO, data.table, MendelianRandomization
- Python packages: matplotlib, numpy (for GTEx validation and figure generation)

## Data Sources
All GWAS summary statistics are publicly available from:
- IEU OpenGWAS: https://gwas.mrcieu.ac.uk/
- eQTLGen: https://eqtlgen.org/
- GTEx v8: https://gtexportal.org/home/datasets
- 1000 Genomes (EUR reference): used for LD clumping

## Script Descriptions

### Phase 2: Univariate MR
- `Phase2_batch_MR.R` — Batch two-sample MR for 8 exposures × 8 outcomes using IVW, MR-Egger, weighted median, weighted mode. Includes heterogeneity, pleiotropy, and MR-PRESSO.

### Phase 3: Multivariable MR
- `Phase3_MVMR.R` — MVMR for 3 groups of correlated exposures (smoking vs alcohol, coffee vs fruit, NO2 vs PM2.5) using MVMR-IVW and MVMR-Egger.

### Phase 5: Three-Step Mediation MR
- `Phase5a_prepare_clump.R` — Step 1: MR of smoking on gene expression (eQTLGen blood cis-eQTL). Prepares LD clumping input.
- `Phase5b_run_clump.sh` — Runs PLINK LD clumping (r² < 0.001, 10 Mb window) using 1000G EUR reference.
- `Phase5c_mediation_MR.R` — Step 2: MR of gene expression on disease. Step 3: Calculates mediation effect and proportion mediated.
- `Phase5d_step3_mediation_effect.R` — Standalone Step 3 script (computes mediation effects from existing Step 1 and Step 2 results).

### Phase 7: Colocalization
- `Phase7_colocalization.R` — Bayesian colocalization (coloc.abf) for 5 exposure-outcome pairs at ±150 kb windows. Includes lead SNP identification, regional association extraction, and PP.H4 calculation.

### Sensitivity Analyses
- `Sensitivity_F_steiger.R` — F-statistics for instrument strength and Steiger directionality test.
- `Sensitivity_leave_one_out.R` — Leave-one-out sensitivity analysis for 14 IVW-significant pairs.
- `Sensitivity_PRESSO_Egger.R` — Compilation of MR-PRESSO, heterogeneity (Cochran's Q), and MR-Egger intercept results.
- `Sensitivity_coloc_50kb.R` — Colocalization sensitivity analysis at ±50 kb window.
- `Sensitivity_GTEx_validation.py` — Multi-tissue eQTL validation using GTEx v8 Lung and Whole Blood data.

## Reproduction Instructions
1. Set OpenGWAS JWT token in `C:/Users/[user]/Documents/.Renviron` as `OPENGWAS_JWT=<token>`
2. Download GTEx v8 eQTL data from https://gtexportal.org/home/datasets
3. Download 1000 Genomes EUR reference panel for LD clumping
4. Run scripts in numerical order (Phase 2 → 3 → 5 → 7 → Sensitivity)
5. Output files are saved to `results/tables/`

## Contact
For questions about the code, please contact the corresponding author.
