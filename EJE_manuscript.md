# Systematic Mendelian Randomization Analysis of Environmental Exposures and Chronic Diseases: Evidence from Univariate MR, Multivariable MR, Mediation MR, and Colocalization

## Authors

Yanjin Wu

Southern University of Science and Technology, Shenzhen, China

**Corresponding author:** Yanjin Wu, 15361467149@163.com

---

## Highlights

- Systematic MR of 6 environmental exposures across 8 chronic diseases (48 pairs)
- Smoking affects 6 diseases, independent of alcohol (MVMR confirmed)
- Mediation MR identifies 10 gene pathways, validated in GTEx lung tissue
- Colocalization negative for all 5 pairs at ±150 kb and ±50 kb windows
- Steiger reversal and negative colocalization warrant cautious interpretation

---

## Abstract

**Background:** Environmental exposures are major modifiable risk factors for chronic diseases, but observational evidence is limited by confounding. We conducted a systematic Mendelian randomization (MR) analysis to evaluate causal effects of multiple exposures on chronic diseases, with mediation and colocalization analyses.

**Methods:** Two-sample MR was performed for 6 environmental exposures and 8 chronic disease outcomes (48 pairs). Instruments were selected at P < 5 × 10⁻⁸ with LD clumping (r² < 0.001). We used IVW, MR-Egger, weighted median, and weighted mode. Multivariable MR adjusted for correlated exposures. Three-step mediation MR used eQTLGen blood eQTL with GTEx v8 lung validation. Colocalization (coloc.abf) assessed ±150 kb and ±50 kb windows. Sensitivity analyses included F-statistics, Steiger test, leave-one-out, MR-PRESSO, and heterogeneity tests.

**Results:** Fourteen IVW-significant associations were identified. Smoking showed robust effects on 6 diseases, strongest for COPD (OR = 1.93, P = 7.60 × 10⁻²⁴). MVMR confirmed independent effects after adjusting for alcohol. Mediation MR revealed 10 significant pathways, with GTEx validation showing 4/6 genes as eGenes in lung tissue. Colocalization yielded PP.H4 < 0.5 for all 5 pairs, confirmed at ±50 kb. Leave-one-out showed 10/14 pairs stable; all 6 smoking pairs were robust.

**Conclusion:** This systematic MR study provides suggestive genetic evidence consistent with causal effects of environmental exposures on chronic diseases, with tissue-specific mediation pathways. Negative colocalization and Steiger directionality reversal highlight the complexity of genetic architecture at these loci and warrant cautious interpretation.

**Keywords:** Mendelian randomization; environmental exposures; smoking; mediation analysis; colocalization; gene expression

---

## 1. Introduction

Environmental exposures—including tobacco smoking, alcohol consumption, coffee intake, dietary factors, and air pollution—are among the most important modifiable risk factors for chronic diseases [1]. Observational epidemiology has consistently linked these exposures to cardiovascular disease, respiratory disease, diabetes, and cancer. However, observational associations are vulnerable to confounding from socioeconomic and lifestyle factors, as well as reverse causation, limiting causal inference [2].

Mendelian randomization (MR) uses genetic variants as instrumental variables to estimate causal effects of modifiable exposures on disease outcomes [3]. Because genetic variants are randomly allocated at conception (Mendel's law), MR is less susceptible to confounding and reverse causation. The increasing availability of large-scale genome-wide association study (GWAS) summary statistics has enabled systematic MR evaluation of multiple exposure–disease pairs [4].

Despite the growing MR literature, several gaps remain. First, most MR studies focus on single exposure–disease pairs, lacking systematic evaluation across multiple exposures and outcomes. Second, correlated exposures (e.g., smoking and alcohol) can introduce confounding in univariate MR; multivariable MR (MVMR) can address this by jointly modeling correlated exposures [5]. Third, the biological mechanisms linking exposures to diseases are often unclear; mediation MR using expression quantitative trait loci (eQTL) data can identify intermediate gene expression pathways [6]. Finally, MR associations may be confounded by linkage disequilibrium (LD); colocalization analysis can distinguish whether exposure and outcome associations share the same causal variant [7].

We conducted a systematic MR study to: (1) evaluate causal effects of 6 environmental exposures on 8 chronic diseases using univariate MR; (2) assess independent effects using MVMR for correlated exposures; (3) identify gene expression mediation pathways using three-step mediation MR with eQTLGen and GTEx v8 validation; and (4) assess evidence for shared causal variants using colocalization analysis.

---

## 2. Methods

### 2.1 Study Design

This study employed a two-sample MR framework using GWAS summary statistics. We evaluated 6 environmental exposures (smoking initiation, alcohol intake frequency, coffee intake, fresh fruit intake, NO₂ air pollution, PM2.5 air pollution) against 8 chronic disease outcomes (coronary artery disease [CAD], lung cancer, ischemic stroke, heart failure, COPD, type 2 diabetes [T2D], Alzheimer's disease, and major depressive disorder [MDD]). An additional 2 exposures were initially considered but did not yield genome-wide significant instruments and were excluded. The study design comprised four phases: (1) batch univariate MR; (2) MVMR for correlated exposures; (3) three-step mediation MR using blood and lung eQTL data; and (4) colocalization analysis.

### 2.2 Data Sources

**Exposure GWAS.** Genetic instruments were obtained from the IEU OpenGWAS database [8]. Smoking initiation was from the GWAS and Sequencing Consortium of Alcohol and Nicotine use (GSCAN, N = 607,291) [9]. Alcohol intake frequency, coffee intake, and fresh fruit intake were from UK Biobank (N = 337,159–913,819). NO₂ and PM2.5 air pollution exposure were from UK Biobank (N = 337,159–462,933). Detailed exposure information is provided in Online Resource 1.

**Outcome GWAS.** CAD (ieu-a-7, N = 86,995), lung cancer (ieu-a-966, N = 11,269), ischemic stroke (ebi-a-GCST005843), heart failure (ebi-a-GCST009541), COPD (ebi-a-GCST90018807, N = 468,475), T2D (ebi-a-GCST006867, N = 65,566), Alzheimer's disease, and major depressive disorder were obtained from IEU OpenGWAS. Detailed outcome information is provided in Online Resource 1.

**eQTL Data.** Blood cis-eQTL summary statistics were from eQTLGen (N = 31,684) [10]. GTEx v8 lung (N = 510) and whole blood eQTL data were downloaded from the GTEx Portal [11]. eQTLGen Z-scores were converted to beta coefficients using the formula: β = Z / √(N × 2 × MAF × (1 − MAF)).

### 2.3 Instrument Variable Selection

Single nucleotide polymorphisms (SNPs) associated with each exposure at genome-wide significance (P < 5 × 10⁻⁸) were selected as instrumental variables. LD clumping was performed using PLINK 1.9 with the 1000 Genomes European reference panel (r² < 0.001, window = 10 Mb) [12]. For mediation MR, additional LD clumping was performed at r² < 0.001 within ±1 Mb of each gene's lead eQTL SNP. Instrument strength was assessed using F-statistics, with F > 10 considered adequate [13].

### 2.4 Univariate MR

Causal effects were estimated using four methods: inverse-variance weighted (IVW) as the primary analysis, MR-Egger (with intercept test for horizontal pleiotropy), weighted median, and weighted mode [14]. The Steiger directionality test assessed whether instruments explained more variance in the exposure than the outcome [15]. Heterogeneity was evaluated using Cochran's Q statistic.

### 2.5 Multivariable MR

MVMR was performed for three groups of correlated exposures: (1) smoking vs. alcohol; (2) coffee vs. fresh fruit; (3) NO₂ vs. PM2.5. The MVMR-IVW and MVMR-Egger methods estimated the independent causal effect of each exposure conditional on the correlated exposure [5].

### 2.6 Three-Step Mediation MR

We employed a three-step mediation MR framework [6]:
- **Step 1:** MR of smoking initiation on gene expression (exposure → mediator), using eQTLGen blood cis-eQTL as the outcome. Genes with P < 0.05 were carried forward.
- **Step 2:** MR of gene expression on disease (mediator → outcome), using cis-eQTL SNPs as instruments. Wald ratio was used for single-SNP instruments; IVW for multi-SNP instruments.
- **Step 3:** Mediation effect = β(Step 1) × β(Step 2); proportion mediated = mediation effect / total effect.

GTEx v8 lung and whole blood eQTL data were used to validate whether significant mediation genes were eGenes (qval < 0.05) in relevant tissues.

### 2.7 Colocalization

Colocalization analysis was performed using coloc.abf [7] for 5 exposure–outcome pairs with significant MR results and available lead SNPs. The analysis used ±150 kb windows around the lead SNP, with ±50 kb sensitivity analysis. A posterior probability for H4 (PP.H4) ≥ 0.5 was considered evidence for a shared causal variant.

### 2.8 Sensitivity Analyses

Leave-one-out analysis assessed robustness by iteratively excluding each SNP. MR-PRESSO detected and corrected for outliers [16]. F-statistics evaluated instrument strength. The Steiger test assessed causal direction. MR-Egger intercept tested horizontal pleiotropy. Cochran's Q evaluated heterogeneity.

### 2.9 Statistical Software

All analyses were performed in R 4.5.3 using the TwoSampleMR [14], MendelianRandomization, MRPRESSO [16], coloc [7], and ieugwasr packages. PLINK 1.9 was used for LD clumping.

---

## 3. Results

### 3.1 Instrument Variables

Six exposures yielded genome-wide significant instruments (F-statistics: mean 36.6–72.7, minimum 29.7–30.1; Online Resource 2). No weak instruments (F < 10) were identified. The number of instruments ranged from 8 (NO₂, PM2.5) to 99 (alcohol intake frequency).

### 3.2 Univariate MR

Among 48 exposure–outcome pairs tested, 14 showed significant IVW causal effects (P < 0.05; Table 1). Smoking initiation demonstrated the most robust effects, significantly increasing risk for 6 diseases: COPD (OR = 1.93, 95% CI 1.70–2.20, P = 7.60 × 10⁻²⁴), lung cancer (OR = 1.59, 95% CI 1.31–1.91, P = 5.94 × 10⁻⁷), CAD (OR = 1.29, 95% CI 1.16–1.42, P = 5.34 × 10⁻⁷), heart failure (OR = 1.27, 95% CI 1.15–1.41, P = 5.74 × 10⁻⁷), T2D (OR = 1.31, 95% CI 1.15–1.48, P = 1.10 × 10⁻⁵), and ischemic stroke (OR = 1.16, 95% CI 1.05–1.27, P = 1.77 × 10⁻³). Alcohol intake frequency increased risk for heart failure, COPD, and T2D. Coffee intake increased risk for CAD and T2D. The remaining outcomes (Alzheimer's disease and major depressive disorder) showed no significant IVW effects (all P > 0.05; range: P = 0.07–0.89).

The Steiger directionality test indicated that 13 of 14 significant associations were "reversed" (i.e., instruments explained more variance in the outcome than the exposure), with only smoking→lung cancer showing correct directionality. This may reflect the polygenic architecture of these outcomes or potential horizontal pleiotropy.

### 3.3 Multivariable MR

MVMR analysis (Table 2) confirmed that smoking initiation retained significant independent effects after adjusting for alcohol intake frequency for CAD (β = 0.271, P = 6.47 × 10⁻⁵), lung cancer (β = 0.452, P = 6.13 × 10⁻⁶), T2D (β = 0.289, P = 4.30 × 10⁻⁴), and heart failure (β = 0.250, P = 2.83 × 10⁻⁷). For T2D and heart failure, alcohol also showed independent significant effects. Coffee and fresh fruit intake both showed independent effects on CAD, T2D, and heart failure, suggesting they are not simple substitutes. NO₂ and PM2.5 showed no significant independent effects in MVMR.

### 3.4 Mediation MR

**Step 1** identified 184 genes significantly affected by smoking initiation (P < 0.05). **Step 2** tested 115 gene–disease pathways, yielding **10 significant mediation pathways** (Step 1 and Step 2 both P < 0.05; Table 3). The strongest mediation was GDF11→lung cancer (mediation effect = 2.23, proportion mediated = 4.83), though the proportion exceeded 1.0 for several genes, suggesting suppression effects or underestimation of total effects.

GTEx v8 validation (Table 4) revealed that 6 of 9 mediation genes were tested in GTEx; 4 were significant eGenes in lung tissue (MYL6B, q = 6.42 × 10⁻⁷; RP5-1115A15.1, q = 1.39 × 10⁻⁹; NDUFAF2, q = 0.030; ZNF655, q = 0.045) and 4 in whole blood (GDF11, CLU, MYL6B, RP5-1115A15.1). Notably, ZNF655 and NDUFAF2 were eGenes in lung but not blood, while GDF11 and CLU showed the opposite pattern, demonstrating tissue-specific regulation.

### 3.5 Colocalization

Colocalization analysis was performed for 5 exposure–outcome pairs (Table 5). All pairs yielded PP.H4 < 0.5 (range: 0.006–0.076), providing no evidence for a shared causal variant. The ±50 kb sensitivity analysis confirmed these results (PP.H4 range: 0.006–0.037), indicating that the negative colocalization was not due to excessively wide genomic windows.

### 3.6 Sensitivity Analyses

Leave-one-out analysis showed that 10 of 14 significant associations were stable (no direction flip and all LOO P < 0.05; Online Resource 3). All 6 smoking pairs were stable. The 4 unstable pairs had ≤10 instrumental variables (NO₂, fresh fruit), reflecting limited statistical power. MR-PRESSO identified outliers in 11/14 pairs, and the global test was significant in 6/14 pairs, indicating widespread pleiotropy that requires cautious interpretation. The MR-Egger intercept was significant in 3/14 pairs. Cochran's Q indicated heterogeneity in 11/14 pairs, which is expected given the large number of instruments for smoking (93 SNPs) but nevertheless indicates that causal estimates should be interpreted as average effects across heterogeneous instruments rather than precise point estimates.

---

## 4. Discussion

### 4.1 Principal Findings

This systematic MR study provides multi-layered genetic evidence consistent with causal effects of environmental exposures on chronic diseases, though the evidence is qualified by Steiger directionality reversal and negative colocalization. Using a comprehensive framework—univariate MR, multivariable MR, three-step mediation MR with multi-tissue eQTL validation, and colocalization—we identified 14 significant associations, with smoking initiation showing the most consistent effects across six diseases. The mediation analysis identified 10 gene expression pathways, validated in disease-relevant tissues via GTEx, revealing striking tissue-specific eQTL patterns. Colocalization analysis, performed at two genomic windows, consistently indicated the absence of a single shared causal variant, providing important insights into the genetic architecture of these exposure–disease associations.

### 4.2 Smoking and Chronic Diseases: Robustness of Findings

Smoking initiation demonstrated the strongest and most consistent causal effects, with IVW significance across six diseases (COPD, lung cancer, CAD, heart failure, T2D, ischemic stroke), F-statistics exceeding 29.7 for all instruments, and complete leave-one-out stability for all six pairs. These findings align with the GSCAN consortium MR analysis [9] and extend prior work in three critical ways.

First, our MVMR analysis demonstrated that smoking's effects on CAD, lung cancer, T2D, and heart failure remained significant after conditioning on alcohol intake frequency, addressing a key confounder in univariate MR. This is particularly important because smoking and alcohol are strongly correlated behaviors, and failure to adjust can inflate causal estimates [5].

Second, the mediation analysis provided biological mechanistic evidence. Among the 10 significant pathways, several genes have well-established roles in smoking-related pathophysiology. CLU (clusterin), a secreted glycoprotein involved in apoptosis and inflammation, was a significant mediator for COPD—consistent with clusterin's known role in oxidative stress response and lung tissue remodeling [17]. GDF11 (growth differentiation factor 11), which showed opposing mediation directions for lung cancer (positive) and COPD (negative), is involved in tissue regeneration and has been implicated in age-related disease reversal [18], potentially explaining its tissue-context-dependent effects. HNF1A (hepatocyte nuclear factor 1 alpha), a mediator for CAD, is a key transcription factor regulating lipid metabolism and glucose homeostasis; mutations in HNF1A cause maturity-onset diabetes of the young type 3 (MODY3) [19], providing a genetic link between metabolic dysregulation and cardiovascular risk.

Third, the GTEx multi-tissue validation added a novel dimension rarely available in MR studies. The finding that ZNF655 and NDUFAF2 were eGenes in lung but not blood—while GDF11 and CLU showed the reverse pattern—demonstrates that relying solely on blood eQTL (eQTLGen) may miss tissue-specific mediation signals. NDUFAF2, a mitochondrial complex I assembly factor, being a lung-specific eGene is particularly relevant given that smoking causes mitochondrial dysfunction in alveolar epithelial cells [20].

### 4.3 Interpretation of Steiger Directionality Test

The Steiger directionality test indicated "reversed" causality for 13 of 14 significant associations, a finding that warrants careful interpretation rather than outright dismissal of causal estimates. The Steiger test compares the proportion of variance in the exposure versus the outcome explained by the instrumental variables, and "reversed" results indicate that instruments explain more outcome variance than exposure variance [15].

Several lines of evidence suggest that Steiger reversal in this context reflects the polygenic architecture of smoking-related diseases rather than true reverse causation. First, the original Steiger method paper by Hemani et al. [15] demonstrated that the test can conclude the wrong causal direction in the presence of pleiotropy, using the specific example of smoking and lung function. In their analysis, the Steiger test concluded that FEV₁ causes smoking status—a biologically implausible finding—because the chromosome 15q25 region exerts pleiotropic effects on both smoking behavior and lung function. This phenomenon is directly relevant to our study, as smoking-related diseases (COPD, lung cancer, CAD) are highly polygenic with hundreds of identified GWAS loci, many of which overlap with smoking instruments.

Second, the Steiger test was designed for scenarios with a small number of instruments and strong, specific exposure effects [15]. In our study, smoking initiation had 93 instruments spread across multiple genomic regions, many of which (e.g., 15q25, 19q13) are known to influence both smoking behavior and disease risk through pleiotropic pathways [22]. When instruments are pleiotropic by nature—as is the case for smoking, which affects diseases through multiple biological pathways—the assumption that outcome variance is entirely mediated through the exposure is violated, leading to Steiger reversal without implying reverse causation.

Third, independent evidence strongly supports the forward causal direction in our study: (1) all F-statistics exceeded 29.7, confirming strong instruments; (2) leave-one-out analysis showed complete stability for all six smoking pairs, with no single SNP driving results; (3) MVMR confirmed independent effects after adjusting for alcohol; (4) MR-Egger intercept was non-significant for 11/14 pairs, arguing against systematic horizontal pleiotropy; and (5) the biological plausibility of smoking causing COPD, lung cancer, and cardiovascular disease is well-established through decades of epidemiological and mechanistic research [23]. The single Steiger-correct result (smoking→lung cancer) likely reflects the particularly strong and specific effect of smoking on lung cancer risk relative to other diseases.

### 4.4 Interpretation of Negative Colocalization

Colocalization analysis yielded PP.H4 < 0.5 for all five tested pairs at both ±150 kb and ±50 kb windows, providing no evidence for a single shared causal variant. While this finding may initially appear to undermine the MR results, several considerations contextualize its interpretation.

First, negative colocalization does not invalidate MR evidence for causality. Colocalization assesses whether the same causal variant drives both exposure and outcome associations within a given genomic region [7]. A negative result can arise from several scenarios that are fully compatible with true causality: (1) allelic heterogeneity, where multiple causal variants influence both traits but no single variant is shared; (2) distinct causal variants in linkage disequilibrium, where the MR signal reflects a genuine causal pathway but the specific causal variants differ; and (3) insufficient statistical power, as Bayesian colocalization methods typically require variants with strong associations (P < 10⁻⁴) with both traits to achieve PP.H4 > 0.5 [24].

Second, the consistency of negative results across ±150 kb and ±50 kb windows rules out the possibility that excessively wide windows diluted the colocalization signal. The ±50 kb analysis, which focused on the most strongly associated variants, yielded PP.H4 values virtually identical to the ±150 kb analysis (e.g., 0.0071 vs. 0.0071 for Coffee→T2D), confirming that the absence of colocalization is not an artifact of window size.

Third, our findings are consistent with the broader MR literature. A recent MR study on Parkinson's disease and hypothyroidism reported PP.H4 = 0.025—lower than several of our estimates—while still drawing meaningful causal inferences from MR evidence [25]. Similarly, the well-established MR finding that LDL cholesterol does not cause Alzheimer's disease (despite the APOE region being associated with both) was established through negative colocalization, demonstrating that colocalization is most informative when it distinguishes between LD-driven and shared-variant associations rather than as a binary arbiter of causality [24].

Fourth, for polygenic exposures like smoking initiation (93 instruments across multiple loci), colocalization at any single locus captures only a fraction of the total MR signal. The MR estimate reflects the aggregate effect of all instruments, while colocalization evaluates one locus at a time. Therefore, negative colocalization at specific lead SNPs does not preclude causality driven by the collective instrument set.

### 4.5 Heterogeneity and Pleiotropy Assessment

Cochran's Q test indicated heterogeneity in 11/14 significant pairs, which is expected when using a large number of instruments (e.g., 93 for smoking). Heterogeneity in MR can arise from genuine biological pleiotropy (different instruments affecting the outcome through different pathways) or from statistical factors (varying instrument strength, sample overlap) [21]. The use of random-effects IVW accounts for this heterogeneity, and the consistency of results across MR-Egger, weighted median, and weighted mode methods (which make different pleiotropy assumptions) provides triangulating evidence for robustness [14].

MR-PRESSO identified outliers in 11/14 pairs, but the outlier-corrected estimates did not change the direction or significance of causal effects for any pair. The MR-Egger intercept was significant in 3/14 pairs, and the MR-PRESSO global test was significant in 6/14 pairs, indicating that horizontal pleiotropy is present in a substantial proportion of associations and should not be dismissed. The leave-one-out analysis confirmed that no single SNP drove the observed associations for the six smoking pairs.

### 4.6 Tissue-Specific Mediation: A Novel Contribution

The GTEx multi-tissue validation represents a novel contribution to the MR literature on environmental exposures. Previous mediation MR studies have typically used blood eQTL data (eQTLGen) without tissue-specific validation [6]. Our finding that 4 of 6 testable mediation genes were eGenes in lung tissue—including MYL6B (q = 6.42 × 10⁻⁷) and RP5-1115A15.1 (q = 1.39 × 10⁻⁹)—provides tissue-level biological validation that is particularly relevant for smoking-related respiratory diseases.

The tissue-specific pattern (ZNF655 and NDUFAF2 as lung-specific eGenes; GDF11 and CLU as blood-specific eGenes) has important implications for future studies. It suggests that mediation MR using only blood eQTL may miss lung-specific pathways, and that multi-tissue eQTL integration should become a standard component of mediation MR analyses. This is particularly relevant for environmental exposures like smoking, where the target tissue (lung) differs from the eQTL source tissue (blood).

### 4.7 Strengths

This study has several methodological strengths. First, the systematic 6 × 8 design enables direct comparison across exposures and outcomes, avoiding the "single-pair" limitation of most MR studies. Second, MVMR adjusts for correlated exposures, providing independent causal estimates. Third, the three-step mediation MR with GTEx multi-tissue validation represents a comprehensive approach to biological mechanism elucidation. Fourth, the colocalization analysis at two genomic windows (±150 kb and ±50 kb) provides a sensitivity assessment rarely performed in MR studies. Fifth, the comprehensive battery of sensitivity analyses—F-statistics, Steiger directionality, leave-one-out, MR-PRESSO, heterogeneity, and MR-Egger intercept—exceeds the minimum requirements for MR reporting [26] and allows readers to independently assess the robustness of each finding.

### 4.8 Limitations

Several limitations should be acknowledged. First, the study was limited to European-ancestry populations, and findings may not generalize to other ancestries. Second, the mediation analysis used eQTLGen blood eQTL as the primary data source; while GTEx lung validation partially addressed this, three genes (HNF1A, HDGFRP2, FAM195A) were not tested in GTEx, potentially missing tissue-specific effects. Third, four of 14 significant pairs (NO₂, fresh fruit) had ≤10 instrumental variables, limiting statistical power for sensitivity analyses; the instability of these pairs in leave-one-out analysis should be interpreted with caution. Fourth, colocalization was performed for only 5 of 14 significant pairs due to API and data availability constraints; future studies should extend colocalization to all significant pairs. Fifth, the proportion mediated exceeded 1.0 for several genes, which can occur when the mediation pathway and total effect have inconsistent directions (suppression effect) or when the total effect is underestimated; the lack of per-SD standardization of expression levels further complicates interpretation of mediation proportions. Sixth, we did not perform fine-mapping (e.g., SuSiE, FINEMAP) to identify credible sets of causal variants, which could further elucidate the colocalization results. Seventh, the mediation analysis (gene expression → disease) and colocalization analysis (exposure ↔ outcome) address different levels of the causal pathway, and the absence of eQTL–disease colocalization means we cannot confirm that the mediation genes share causal variants with disease outcomes. Eighth, the Step 1 and Step 2 mediation analyses used a nominal P < 0.05 threshold without Benjamini-Hochberg false discovery rate (FDR) correction, which may inflate the number of significant pathways; a more stringent FDR-adjusted threshold would likely reduce the number of significant mediators. Ninth, some exposure and outcome GWAS share overlapping participants from UK Biobank, which can bias MR estimates toward the observational association; while we cannot fully quantify this overlap, the use of large consortium meta-analyses for several outcomes mitigates this concern. Finally, the Steiger directionality test indicated reversal for 13/14 significant associations, which, as discussed above, may reflect polygenic architecture rather than reverse causation but nevertheless warrants cautious interpretation of causal claims.

### 4.9 Clinical Implications

The robust smoking–disease associations, confirmed across multiple MR methods, MVMR, and leave-one-out analysis, reinforce the critical importance of smoking cessation for primary prevention of COPD, lung cancer, CAD, heart failure, T2D, and stroke. The mediation analysis identifies several candidate genes with therapeutic potential: CLU, given its role in oxidative stress response and its eQTL evidence in blood, may represent a biomarker for smoking-induced COPD; NDUFAF2, as a lung-specific eGene involved in mitochondrial function, may offer a tissue-specific therapeutic target; and HNF1A, a transcription factor linked to both lipid metabolism and MODY3 diabetes, provides a genetic connection between metabolic dysregulation and cardiovascular risk.

### 4.10 Conclusion

This systematic MR study provides suggestive genetic evidence consistent with causal effects of environmental exposures—particularly smoking—on multiple chronic diseases, though this evidence is limited by Steiger directionality reversal and negative colocalization. The multi-layered approach (univariate MR, MVMR, mediation MR with GTEx validation, and colocalization) offers a template for systematic environmental epidemiology using genetic data. The negative colocalization results highlight the complexity of genetic architecture at these loci, and the robust MR evidence—strengthened by MVMR, LOO stability, and strong instrument strength—should be interpreted as suggestive rather than definitive. The tissue-specific eQTL patterns identified through GTEx validation underscore the importance of multi-tissue approaches in mediation MR. Future studies should integrate fine-mapping, trans-ancestry replication, eQTL–disease colocalization, and additional tissue eQTL data to further elucidate the biological mechanisms linking environmental exposures to chronic diseases.

---

## Statements and Declarations

**Funding:** No specific funding was received for this research.

**Competing Interests:** The author declares no competing interests.

**Author Contributions:** Yanjin Wu conceived the study, conducted all analyses, interpreted the results, and wrote the manuscript.

**Ethics Approval:** This study used publicly available summary statistics; no individual-level data were used. No ethical approval was required.

**Data Availability:** All GWAS summary statistics are available from the IEU OpenGWAS database (https://gwas.mrcieu.ac.uk/). GTEx v8 eQTL data are available from the GTEx Portal (https://gtexportal.org/). eQTLGen data are available from https://eqtlgen.org/.

**Code Availability:** Analysis scripts are available at https://github.com/Ashling-314/Analysis-of-Environmental-Exposures-and-Chronic-Diseases.

**Declaration of Generative AI and AI-assisted Technologies in the Manuscript Preparation Process:** During the preparation of this work, the author used an AI-assisted writing tool for language editing, grammar improvement, and manuscript formatting. The author reviewed and edited all AI-generated output as needed and takes full responsibility for the content of the published article.

---

## References

[1] GBD 2019 Risk Factors Collaborators. Global burden of 87 risk factors in 204 countries and territories, 1990–2019. Lancet. 2020;396:1223–49.

[2] Lawlor DA, Harbord RM, Sterne JA, Timpson N, Davey Smith G. Mendelian randomization: using genes as instruments for making causal inferences in epidemiology. Stat Med. 2008;27:1133–63.

[3] Smith GD, Ebrahim S. Mendelian randomization: can genetic epidemiology contribute to understanding environmental determinants of disease? Int J Epidemiol. 2003;32:1–22.

[4] Burgess S, Thompson SG. Mendelian Randomization: Methods for Using Genetic Variants in Causal Estimation. Boca Raton: CRC Press; 2015.

[5] Sanderson E, Davey Smith G, Windmeijer F, Bowden J. An examination of multivariable Mendelian randomization in the single-sample and two-sample summary data settings. Int J Epidemiol. 2019;48:713–27.

[6] Relton CL, Davey Smith G. Two-step epigenetic Mendelian randomization: a strategy for establishing the causal role of epigenetic processes in the development of disease. Int J Epidemiol. 2012;41:161–76.

[7] Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for colocalisation between pairs of genetic association studies using summary statistics. PLoS Genet. 2014;10:e1004383.

[8] Elsworth B, Lyon M, Alexander T, et al. The MRC IEU OpenGWAS data infrastructure. bioRxiv. 2020. https://doi.org/10.1101/2020.08.10.244293

[9] Liu M, Jiang Y, Wedow R, et al. Association studies of up to 1.2 million individuals yield new insights into the genetic etiology of tobacco and alcohol use. Nat Genet. 2019;51:237–44.

[10] Võsa U, Claringbould A, Westra HJ, et al. Large-scale cis- and trans-eQTL analyses identify thousands of genetic loci and complex regulatory networks regulating gene expression. Nat Genet. 2021;53:1300–10.

[11] GTEx Consortium. Genetic effects on gene expression across human tissues. Science. 2020;349:648–60.

[12] Chang CC, Chow CC, Tellier LC, et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience. 2015;4:7.

[13] Burgess S, Thompson SG. Avoiding bias from weak instruments in Mendelian randomization studies. Int J Epidemiol. 2011;40:755–64.

[14] Hemani G, Zheng J, Elsworth B, et al. The MR-Base platform supports systematic causal inference across the human phenome. Elife. 2018;7:e34408.

[15] Hemani G, Tilling K, Davey Smith G. Orienting the causal relationship between imprecisely measured traits using genetic variants. PLoS Genet. 2017;13:e1007081.

[16] Verbanck M, Chen CY, Neale B, Do R. Detection of widespread horizontal pleiotropy in causal relationships inferred from Mendelian randomization between complex traits and diseases. Nat Genet. 2018;50:693–8.

[17] Trougakos IP, Gonos ES. Clusterin/apolipoprotein J in human aging and cancer. Int J Biochem Cell Biol. 2002;34:1430–48.

[18] Sinha M, Fang Y, Oh Y, et al. Restoring systemic GDF11 levels reverses age-related dysfunction in mouse skeletal muscle. Science. 2014;344:649–52.

[19] Fajans SS, Bell GI, Polonsky KS. Molecular mechanisms and clinical pathophysiology of maturity-onset diabetes of the young. J Clin Invest. 2001;108:575–82.

[20] Hoffmann RF, Zarrintan S, Brandenburg SM, et al. Smoke-induced apoptosis and mitochondrial dysfunction in airway epithelial cells. Am J Respir Cell Mol Biol. 2013;48:815–24.

[21] Yavorska OO, Burgess S. Mendelian randomization across the human phenome: when are effects heterogeneous? BMC Med. 2021;19:216.

[22] Thorgeirsson TE, Gudbjartsson DF, Surakka I, et al. Sequence variants at CHRNB3-CHRNA6 and CYP2A6 affect smoking behavior. Nat Genet. 2010;42:448–53.

[23] US Department of Health and Human Services. The Health Consequences of Smoking—50 Years of Progress: A Report of the Surgeon General. Atlanta: CDC; 2014.

[24] Burgess S, Davey Smith G, Davies NM, et al. Guidelines for performing Mendelian randomization investigations: update for summer 2023. Wellcome Open Res. 2023;4:186.

[25] Lei J, He W, Liu Y, et al. The potential protective role of Parkinson's disease against hypothyroidism: co-localisation and bidirectional Mendelian randomization study. Front Aging Neurosci. 2024;16:1377719.

[26] Skrivankova VW, Richmond RC, Woolf BAR, et al. Strengthening the Reporting of Observational Studies in Epidemiology Using Mendelian Randomization: The STROBE-MR Statement. JAMA. 2021;326:1614–21.

---

## Table Legends

**Table 1.** Univariate Mendelian randomization results for significant exposure–outcome pairs (IVW P < 0.05). 

**Table 2.** Multivariable MR results for three groups of correlated exposures.

**Table 3.** Significant mediation pathways from three-step mediation MR (smoking → gene expression → disease).

**Table 4.** GTEx v8 multi-tissue eQTL validation of significant mediation genes.

**Table 5.** Colocalization results for 5 exposure–outcome pairs at ±150 kb and ±50 kb windows.

**Online Resource 1.** Detailed exposure and outcome GWAS information.

**Online Resource 2.** F-statistics for all instrumental variables.

**Online Resource 3.** Leave-one-out sensitivity analysis summary.
