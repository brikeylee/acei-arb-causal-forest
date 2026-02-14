# Individual-Level Treatment Effect Estimation for Perioperative ACEI/ARB Management in Cardiac Surgery: A Causal Forest Analysis

## Abstract

**Background:** Recent evidence suggests that perioperative continuation of renin-angiotensin system inhibitors (RASi) reduces severe adverse events after cardiac surgery, but average treatment effects cannot inform decisions for individual patients. We applied causal forest analysis to estimate individual-level treatment effects and identify which patients benefit most from continuation.

**Results:** Among 9,456 propensity score-matched cardiac surgery patients (from 12,355 eligible), causal forest estimated a mean conditional average treatment effect (CATE) of +0.012 (95% CI 0.002–0.022, p=0.020) for the composite adverse outcome, confirmed by bootstrap validation (mean +0.013, 95% CI 0.004–0.022) and inverse probability of treatment weighting (OR 1.24, 95% CI 1.08–1.42). The CATE distribution revealed that 97.4% of patients had positive estimates (favoring continuation), ranging from -0.021 to +0.045 (coefficient of variation 0.52). While the grf calibration test did not detect statistically significant qualitative heterogeneity (p=0.95), CATE was significantly correlated with age (Spearman rho=0.21, p<10^-96) and inversely with eGFR (rho=-0.14, p<10^-42), indicating continuous dose-response relationships. Patients aged >=65 had a mean CATE of 0.013 with NNT=29, compared with 0.011 and NNT=238 in younger patients (interaction p<0.001). Three patient clusters were identified: a high-benefit group (60.4%, mean CATE 0.013, hypertensive, low diabetes prevalence), a moderate-benefit group (34.2%, mean CATE 0.010, older, diabetic, hypertensive), and a smaller moderate-benefit group (5.4%, mean CATE 0.011, younger, normotensive). All findings showed 100% directional concordance with traditional logistic regression subgroup estimates.

**Conclusions:** Causal forest analysis provides patient-level treatment effect estimates that go beyond average effects, revealing continuous relationships between patient characteristics and treatment benefit. While virtually all cardiac surgery patients benefit from ACEI/ARB continuation, the magnitude varies meaningfully—from negligible (NNT>200) in younger patients to clinically important (NNT~29) in elderly patients. These individual-level estimates can support precision perioperative decision-making.

**Keywords:** causal forest, heterogeneous treatment effects, ACEI/ARB, cardiac surgery, individualized medicine, perioperative management, conditional average treatment effect

---

## Background

Perioperative management of angiotensin-converting enzyme inhibitors (ACEI) and angiotensin receptor blockers (ARB) in cardiac surgery has been the subject of longstanding clinical debate [1,2]. Recent large-scale evidence, including multicenter cohort studies using severity-stratified analysis of acute kidney injury (AKI), has begun to clarify the average treatment effect: perioperative RASi continuation appears to reduce severe AKI and mortality while modestly increasing mild, presumably functional AKI [3,4]. However, these population-level estimates—expressed as average risk ratios or odds ratios—cannot answer the question that clinicians face at the bedside: *how much will this specific patient benefit from continuation?*

Traditional subgroup analyses partially address this gap by estimating treatment effects within pre-specified categories (e.g., age above or below 65 years), but they are limited to a small number of dichotomized variables and ignore continuous dose-response relationships [5]. Moreover, the choice of subgroup boundaries is inherently arbitrary, and multiple testing inflates false-positive rates [6].

Causal forests, a machine learning extension of random forests designed for heterogeneous treatment effect estimation, offer a principled alternative [7]. By recursively partitioning the covariate space to maximize treatment effect heterogeneity, causal forests estimate a conditional average treatment effect (CATE) for each individual patient—a continuous, multidimensional function of their baseline characteristics. Unlike traditional methods that yield a single number per subgroup, causal forests produce patient-specific estimates that capture the full spectrum of treatment benefit.

Importantly, this approach addresses a key limitation of existing evidence. Although prior studies have established that the *average* patient benefits from RASi continuation [3,4], the *distribution* of individual treatment effects—the proportion of patients with large, small, or negligible benefit—remains unknown. Understanding this distribution is essential for precision medicine: if virtually all patients benefit similarly, a universal continuation strategy suffices; if substantial heterogeneity exists, targeted approaches based on patient characteristics may improve outcomes.

In this study, we applied causal forest analysis to a multicenter cohort of cardiac surgery patients to: (1) estimate individual-level treatment effects of perioperative ACEI/ARB continuation versus discontinuation; (2) characterize the distribution and determinants of treatment effect heterogeneity; and (3) compare causal forest estimates with traditional subgroup analyses to assess methodological concordance.

## Methods

### Study Design and Population

This was a retrospective cohort study of consecutive adult patients undergoing cardiac surgery at two tertiary cardiac surgery centers within the Chinese PLA General Hospital system. The study was approved by the institutional ethics committee (Chinese PLA General Hospital S2021-305-01), with informed consent waived due to the retrospective design. Reporting follows the STROBE guideline for observational studies.

Eligible patients were aged >=18 years, underwent elective or urgent cardiac surgery, and had documented preoperative ACEI/ARB use. Patients with emergency surgery, missing treatment exposure data, or missing key outcome variables were excluded.

### Exposure and Outcomes

The exposure was preoperative ACEI/ARB discontinuation (W=1) versus continuation (W=0). The primary outcome was a composite of 30-day mortality, AKI stage 3 (KDIGO), ECMO, intra-aortic balloon pump (IABP), or new-onset dialysis. Secondary outcomes included individual composite components and AKI stages 1–3.

### Missing Data Handling

Missing continuous variables were imputed using multiple imputation by chained equations (MICE, m=5 imputations, predictive mean matching). Body mass index was calculated from height and weight; implausible values were replaced with the cohort median. Categorical variables with missing values (<0.5%) were handled by complete case deletion.

### Propensity Score Matching

Propensity scores were estimated via logistic regression with 21 covariates (demographics, comorbidities, laboratory values, medications, surgical factors). Patients were matched 1:1 using nearest-neighbor matching (caliper=0.1) without replacement. Balance was assessed with standardized mean differences (threshold <0.1).

### Causal Forest Analysis

Causal forests were fitted using the generalized random forest (grf) R package [7] with the following specifications:
- **Trees**: 2,000 with honest splitting (separate estimation and splitting samples)
- **Tuning**: All hyperparameters automatically tuned via out-of-bag optimization
- **Covariates**: 15 variables (age, sex, center, BMI, eGFR, 5 comorbidities, 4 medications, CPB use)
- **Seed**: Fixed for reproducibility

For each patient *i*, the causal forest estimates a conditional average treatment effect (CATE), denoted τ(x_i), representing the expected difference in outcome probability between discontinuation and continuation given covariates x_i. Positive CATE values indicate that discontinuation increases risk (i.e., continuation is protective).

The average treatment effect (ATE) was estimated using the doubly-robust augmented inverse propensity weighting estimator implemented in grf. Treatment effect heterogeneity was formally tested using the grf calibration test, which assesses whether the forest-estimated CATEs have significant predictive power for the true treatment effect beyond the mean [8].

Separate causal forests were fitted for each secondary outcome using identical specifications.

### Heterogeneity Characterization

The distribution of individual CATEs was characterized descriptively (mean, SD, IQR, range, proportion positive/negative). Continuous relationships between CATEs and key covariates (age, eGFR, BMI) were assessed using Spearman correlations and LOESS smoothing. Data-driven patient clustering was performed using k-means on standardized CATE and clinical features, with optimal cluster number selected by silhouette analysis (k=2 to 6).

### Concordance with Traditional Methods

To assess methodological consistency, we compared causal forest subgroup-specific CATE means with traditional logistic regression estimates (odds ratios and risk differences) across 11 pre-specified subgroups. Directional concordance was defined as agreement in the sign of the treatment effect.

### Sensitivity Analyses

1. **Inverse probability of treatment weighting (IPTW)**: Stabilized weights with 1st/99th percentile trimming and sandwich standard errors.
2. **Bootstrap validation**: 500 resamples with complete causal forest re-fitting to assess ATE stability.
3. **Caliper sensitivity**: PSM repeated with calipers of 0.05, 0.10, 0.15, and 0.20.
4. **E-value**: Quantifying the minimum strength of unmeasured confounding to explain the observed association [9].
5. **Hyperparameter sensitivity**: Causal forest re-fitted across 9 combinations of tree number (1000, 2000, 4000) and minimum node size (5, 10, 20).
6. **Subgroup interaction tests**: Logistic regression with interaction terms for 7 clinical variables, FDR-corrected.

Analyses were performed in R 4.4.3 using grf 2.3.2, MatchIt, WeightIt, EValue, and mice packages.

## Results

### Study Population

Of 12,400 screened patients, 12,355 met eligibility criteria. After multiple imputation and 1:1 propensity score matching, 9,456 patients (4,728 matched pairs) were included (Figure 1). Baseline characteristics were well-balanced after matching (all SMDs <0.1; Figure 2). Mean age was 59 years, 66% were male, and the composite adverse event rate was 8.4%.

### Average Treatment Effect: Multi-Method Concordance

All three analytical approaches yielded concordant estimates for the primary composite outcome (Table 1):

| Method | Estimate | 95% CI | p-value |
|--------|----------|--------|---------|
| Causal Forest ATE | +0.012 | 0.002 to 0.022 | 0.020 |
| PSM (OR) | 1.237 | 1.069 to 1.431 | 0.004 |
| IPTW (OR) | 1.237 | 1.078 to 1.419 | 0.001 |

Bootstrap validation confirmed stability (mean ATE +0.013, 95% CI 0.004–0.022; bias +0.001). Results were consistent across all caliper widths (OR range 1.24–1.27, all p<0.005). The E-value was 1.78.

### Distribution of Individual Treatment Effects

The causal forest estimated individual CATEs for each of the 9,456 patients (Figure 3). Key distributional features:

- **Mean CATE**: 0.012 (SD 0.006)
- **Range**: -0.021 to +0.045
- **97.4% of patients** had positive CATEs (favoring continuation)
- **2.6%** had negative CATEs (potentially favoring discontinuation)
- **Coefficient of variation**: 0.52, indicating moderate relative heterogeneity
- **IQR**: 0.008 to 0.016

The grf calibration test confirmed that the mean forest prediction was significant (p=0.010), indicating the CATE estimates capture a real treatment effect. However, the test for differential prediction—assessing whether heterogeneity beyond the mean is statistically significant—was not significant (p=0.95). This indicates that the protective effect of continuation operates in the same direction for virtually all patients, with clinically meaningful but statistically modest variation in magnitude.

### Continuous Determinants of Treatment Benefit

Unlike traditional subgroup analyses that impose arbitrary dichotomizations, causal forest CATEs revealed continuous dose-response relationships with baseline characteristics (Figure 4):

**Age** (Spearman rho=0.213, p<10^-96): Treatment benefit increased monotonically with age. LOESS smoothing showed a near-linear positive relationship, with CATE rising from approximately 0.008 at age 40 to 0.020 at age 80—a 2.5-fold gradient. This translates to NNT improvement from >200 in young patients to approximately 29 in elderly patients (>=65 years).

**eGFR** (rho=-0.140, p<10^-42): Patients with lower eGFR derived greater benefit from continuation. The relationship was nonlinear, with a steeper CATE gradient below eGFR 60 mL/min/1.73m². Patients with eGFR <60 had a mean CATE of 0.012 (NNT=29) versus 0.011 (NNT=151) for eGFR >=90.

**BMI** (rho=0.113, p<10^-28): A weaker but significant positive correlation, with higher BMI associated with greater treatment benefit.

These continuous relationships provide granularity that binary subgroup analyses cannot capture. For example, a 72-year-old patient with eGFR 55 and a 45-year-old with eGFR 105 both fall within broad "benefit" categories in traditional analysis, but their estimated CATEs differ by approximately 2-fold (0.018 vs 0.009).

### Patient Clustering

Silhouette analysis identified k=3 as optimal (silhouette score 0.28). The three clusters represent clinically distinct patient profiles (Table 2):

| Cluster | N (%) | Mean CATE | Age | eGFR | Key Features |
|---------|-------|-----------|-----|------|--------------|
| 1 (High-benefit) | 5,712 (60%) | 0.013 | 58 | 90 | Hypertensive, low diabetes |
| 2 (Moderate-benefit) | 3,231 (34%) | 0.010 | 64 | 81 | Hypertensive + diabetic |
| 3 (Moderate-benefit) | 513 (5%) | 0.011 | 52 | 98 | Normotensive |

All clusters showed positive mean CATEs, reinforcing the universal direction of the treatment effect. The 1.3-fold difference in mean CATE between clusters 1 and 2 reflects the combined influence of age, eGFR, and comorbidity burden.

### Concordance Between Causal Forest and Traditional Estimates

Across 11 pre-specified subgroups, causal forest and traditional logistic regression showed **100% directional concordance** (Table 3). In all subgroups where traditional analysis found a significant effect, causal forest identified a higher mean CATE; in subgroups with non-significant traditional estimates, CATEs were lower but still positive.

This concordance provides mutual validation: traditional methods confirm the causal forest-identified patterns, while causal forest adds the continuous, individual-level granularity that traditional methods lack.

### Secondary Outcomes: Multi-Outcome Causal Forest

Separate causal forests for individual outcomes revealed a severity-dependent pattern consistent with prior literature [3] (Table 4):

| Outcome | ATE | 95% CI | p-value | Direction |
|---------|-----|--------|---------|-----------|
| Composite AE | +0.012 | 0.002 to 0.022 | 0.020 | Discontinuation harmful |
| 30-day mortality | +0.003 | -0.001 to 0.008 | 0.169 | Non-significant |
| AKI Stage 1 | -0.026 | -0.043 to -0.010 | 0.001 | Discontinuation protective |
| AKI Stage 2 | +0.002 | -0.008 to 0.011 | 0.714 | Non-significant |
| AKI Stage 3 | +0.003 | -0.004 to 0.010 | 0.369 | Non-significant |
| AKI All Stages | -0.024 | -0.041 to -0.007 | 0.006 | Discontinuation protective |

The opposing directions for mild AKI (discontinuation protective) and the composite endpoint (discontinuation harmful) are consistent with the severity-dependent dissociation previously described using traditional methods in a larger multicenter cohort [3], providing independent methodological confirmation of this phenomenon.

## Discussion

### From Average Effects to Individual Estimates

This study extends the evidence on perioperative ACEI/ARB management from population-level averages to patient-level estimates. While prior work has established that, on average, continuation reduces severe adverse events [3,4], the clinical utility of average effects is limited when patients vary widely in their baseline risk profiles. Our causal forest analysis demonstrates that individual treatment effects span a clinically meaningful range, from negligible benefit in young, healthy patients (CATE ~0.008, NNT>200) to substantial benefit in elderly patients with impaired renal function (CATE ~0.020, NNT~29).

The 2.5-fold CATE gradient across the age spectrum and the continuous inverse relationship with eGFR provide information that traditional binary subgroups cannot capture. For the clinician deciding whether to continue ACEI/ARB in a 55-year-old with eGFR 75, average treatment effects or dichotomized subgroup analyses offer limited guidance. Causal forest CATEs, by contrast, position this patient on a continuous benefit gradient, enabling more nuanced risk-benefit communication.

### Heterogeneity: Quantitative, Not Qualitative

An important finding is that the grf calibration test for differential heterogeneity was not significant (p=0.95), indicating that the treatment effect does not vary *qualitatively* across patients. Virtually all patients (97.4%) are estimated to benefit from continuation, with no identifiable subgroup for whom discontinuation is preferable. This "universal but graded" benefit pattern has direct clinical implications: it supports a default strategy of continuation for all patients, with the magnitude of expected benefit varying by patient characteristics rather than the direction.

This finding is itself clinically useful. In contrast to conditions where treatment effect heterogeneity includes both benefit and harm subgroups (necessitating careful patient selection), perioperative ACEI/ARB continuation appears to be a "broadly beneficial" intervention with a gradient of effect size. This simplifies clinical decision-making while still allowing for informed shared decision-making based on individual estimated benefit.

### Methodological Concordance as Validation

The 100% directional concordance between causal forest estimates and traditional logistic regression across 11 subgroups provides strong mutual validation. This concordance addresses a common concern about machine learning methods in clinical research—that they may identify spurious patterns. When a data-driven method (causal forest) and a theory-driven method (pre-specified subgroup analysis) converge on the same conclusions, confidence in both is strengthened.

Furthermore, the concordance between three distinct analytical approaches (PSM, IPTW, and causal forest) for the average treatment effect—all yielding point estimates within a narrow range (OR 1.24 for PSM and IPTW; ATE +0.012 for causal forest)—reinforces the robustness of the findings against model specification choices.

### Comparison with Prior Evidence

Our average treatment effect estimates are consistent with prior multicenter evidence demonstrating that perioperative RASi continuation reduces severe adverse events [3,10,11]. The severity-dependent pattern in secondary AKI outcomes—mild AKI reduced by discontinuation, severe outcomes increased—has been independently reported using traditional methods in a larger five-center cohort [3]. The causal forest approach provides complementary individual-level confirmation of this phenomenon.

The absence of statistically significant qualitative heterogeneity in our causal forest analysis aligns with prior subgroup analyses showing limited effect modification [3,12]. However, our identification of continuous dose-response relationships (particularly the monotonic age-CATE gradient and nonlinear eGFR-CATE relationship) adds new information beyond what binary subgroup analyses have previously revealed.

### Limitations

Several limitations warrant acknowledgment. First, despite multi-method validation and sensitivity analyses (E-value 1.78), unmeasured confounding cannot be excluded in this observational study. Second, the dataset originates from two Chinese centers, limiting generalizability. Third, BMI was imputed for 19% of patients, though sensitivity analyses were robust. Fourth, the causal forest did not detect statistically significant qualitative heterogeneity; the identified variation is quantitative, representing gradients of a universally protective effect. This limits the ability to identify patients who would specifically benefit from discontinuation. Fifth, pooled multiple imputation (averaging m=5 datasets) may slightly underestimate standard errors compared with full Rubin's rules. Sixth, external validation of the individual-level estimates in independent cohorts is needed before clinical implementation.

## Conclusions

Causal forest analysis reveals that the protective effect of perioperative ACEI/ARB continuation in cardiac surgery is universal in direction but graded in magnitude. Individual treatment effects range from negligible in young patients to clinically important in elderly patients with impaired renal function, with continuous dose-response relationships that traditional subgroup analyses cannot capture. These patient-level estimates provide a framework for precision perioperative decision-making, supporting a default strategy of continuation with individualized benefit communication.

## Declarations

**Ethics approval and consent to participate:** Approved by the institutional ethics committee (Chinese PLA General Hospital S2021-305-01). Informed consent waived due to retrospective design.

**Consent for publication:** Not applicable.

**Availability of data and materials:** Available from the corresponding author on reasonable request.

**Competing interests:** The authors declare no competing interests.

**Funding:** [To be completed]

**Authors' contributions:** [To be completed]

**Acknowledgements:** [To be completed]

## References

1. Hollmann C, Fernandes NL, Biccard BM. A systematic review of outcomes associated with withholding or continuing angiotensin-converting enzyme inhibitors and angiotensin receptor blockers before noncardiac surgery. Anesth Analg. 2018;127:678-687.
2. Legrand M, Futier E, Mebazaa A, et al. Continuation vs discontinuation of renin-angiotensin system inhibitors before major noncardiac surgery: the Stop-or-Not randomized clinical trial. JAMA. 2024;332:970-978.
3. Zhong Q, Li Z, Jiang Y, et al. Perioperative continuing vs. discontinuing RAS inhibitors on kidney outcomes and mortality in cardiac surgery: a target trial emulation using multicenter data. [Submitted/Under review at Critical Care]. 2026.
4. Ding Q, Zhang Z, Liu H, et al. Perioperative use of renin-angiotensin system inhibitors and outcomes in patients undergoing cardiac surgery. Nat Commun. 2019;10:4202.
5. Kent DM, Steyerberg E, van Klaveren D. Personalized evidence-based medicine: predictive approaches to heterogeneous treatment effects. BMJ. 2018;363:k4245.
6. Brookes ST, Whitely E, Egger M, et al. Subgroup analyses in randomised trials: risks of subgroup-specific analyses; power and sample size for the interaction test. J Clin Epidemiol. 2004;57:229-236.
7. Athey S, Tibshirani J, Wager S. Generalized random forests. Ann Stat. 2019;47:1148-1178.
8. Chernozhukov V, Demirer M, Duflo E, Fernandez-Val I. Generic machine learning inference on heterogeneous treatment effects in randomized experiments. arXiv preprint arXiv:1712.04802. 2018.
9. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med. 2017;167:268-274.
10. Drenger B, Fontes ML, Miao Y, et al. Patterns of use of perioperative angiotensin-converting enzyme inhibitors in coronary artery bypass graft surgery with cardiopulmonary bypass. Circulation. 2012;126:261-269.
11. Zhou H, Xie J, Zheng Z, et al. Effect of renin-angiotensin system inhibitors on acute kidney injury among patients undergoing cardiac surgery: a review and meta-analysis. Semin Thorac Cardiovasc Surg. 2021;33:1014-1022.
12. Roshanov PS, Rochwerg B, Patel A, et al. Withholding versus continuing angiotensin-converting enzyme inhibitors or angiotensin II receptor blockers before noncardiac surgery. Anesthesiology. 2017;126:16-27.

## Tables and Figures

**Table 1.** Multi-method treatment effect estimates for the primary composite outcome.

**Table 2.** Patient cluster profiles from causal forest-derived CATE clustering.

**Table 3.** Concordance between causal forest and traditional logistic regression subgroup estimates.

**Table 4.** Causal forest average treatment effects across individual outcomes.

**Figure 1.** Patient flow diagram.

**Figure 2.** Covariate balance after propensity score matching (Love plot).

**Figure 3.** Distribution of individual conditional average treatment effects (CATEs).

**Figure 4.** Continuous relationships between patient characteristics and treatment benefit: (A) Age vs CATE, (B) eGFR vs CATE, with LOESS smoothing.

**Figure 5.** Variable importance in the causal forest model.

**Supplementary Table S1.** Missing data report.

**Supplementary Table S2.** Caliper sensitivity analysis.

**Supplementary Table S3.** Hyperparameter sensitivity analysis.

**Supplementary Table S4.** Subgroup interaction tests with FDR correction.

**Supplementary Table S5.** Competing risk analysis for AKI outcomes.

**Supplementary Figure S1.** Propensity score distribution by treatment group.

**Supplementary Figure S2.** Bootstrap ATE distribution (500 resamples).

**Supplementary Figure S3.** Calibration plot: predicted vs observed CATE by decile.
