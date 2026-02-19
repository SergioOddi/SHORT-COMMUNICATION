# Discussion — Notes and Critical Points

Questo file raccoglie osservazioni critiche e punti da sviluppare nella Discussion del manoscritto.

---

## D.1 Female Tg2576 mice are cognitively impaired but unresponsive to PEA

A key finding is that the lack of PEA efficacy in females cannot be attributed to an absence of disease. Vehicle-treated female Tg2576 mice displayed clear cognitive impairment relative to wild-type littermates (NORT: 53.8% vs. 66.1%), and elevated 3-NT levels indicative of active nitrosative stress. The amyloid-driven pathology is therefore present in both sexes at 12 months, yet PEA selectively rescues males.

This dissociation is critical for interpretation: it rules out the trivial explanation that "females don't need treatment because they are not impaired" and instead points to a sex-specific biological factor that limits PEA therapeutic activity. The consistent attenuation of effect sizes across three independent endpoints (NOR, spine density, 3-NT) strengthens this conclusion, as it is unlikely that three unrelated ceiling effects would produce the same pattern.

## D.2 Sexual dimorphism in baseline measures

### 3-Nitrotyrosine
Male Tg2576 mice show approximately twice the 3-NT levels of females (19.1 vs. 9.7), indicating a pronounced sexual dimorphism in nitrosative stress burden. This likely reflects known sex differences in oxidative/nitrosative vulnerability, potentially modulated by estrogen-mediated antioxidant pathways. However, the female 3-NT levels remain elevated relative to wild-type controls, confirming ongoing pathological stress.

Two interpretive considerations:
- The lower baseline in females reduces the absolute "room for improvement," potentially contributing to the smaller observed effect. However, the relative reduction by PEA is also smaller in females (~21% vs. ~30% in males), suggesting that baseline differences alone do not fully explain the attenuated response.
- The near-significant Sex × Treatment interaction (p = 0.059) in the ANOVA supports the interpretation that PEA efficacy genuinely differs between sexes, rather than reflecting a simple scaling effect.

### Dendritic spine density
Female Tg2576 mice exhibited substantially higher baseline spine densities than males in both CA1 compartments (apical: 2.21 vs. 1.42 spines/μm; basal: 2.11 vs. 1.21). This difference is not disease-specific, as it was also observed in wild-type littermates (apical: 2.30 vs. 1.66; basal: 2.26 vs. 1.47), indicating an intrinsic sexual dimorphism in hippocampal dendritic morphology.

This raises the ceiling effect hypothesis: if females retain higher synaptic density despite the transgene, the margin for PEA-mediated improvement may be structurally limited. Under this interpretation, the lack of PEA effect on spines in females would partly reflect reduced room for improvement.

However, this explanation is specific to the spine endpoint and cannot account for the attenuated PEA effect on cognition (NORT), where females do show clear impairment. The dissociation between cognitive deficit (present) and spine preservation (partial) in female Tg2576 also suggests that spine density alone may not fully capture the synaptic dysfunction underlying memory impairment — functional synaptic alterations (e.g., LTP impairment, receptor changes) may contribute independently.

## D.3 Compartment-specific PEA effect on spines (apical vs. basal)

PEA significantly increased basal but not apical spine density in males. This compartment dissociation was not reported by Tortolani et al. (2025), who found effects on both compartments. Possible explanations:

1. **Statistical approach.** The original study used ANOVA with post-hoc tests across all 8 groups (WT + Tg × Sex × Treatment), while our planned contrasts focus specifically on within-sex Tg comparisons. Different analytic frameworks may yield different significance thresholds.

2. **Biological heterogeneity.** Apical and basal dendrites receive distinct afferent inputs (Schaffer collaterals vs. commissural/associational fibers) and may be differentially vulnerable to Aβ-mediated toxicity and differentially responsive to PEA-mediated neuroprotection.

3. **Sample variability.** The apical male PEA group shows notably high variability (IQR = 0.78 vs. 0.40 in vehicle), which reduces statistical power. This variability may reflect biological heterogeneity in treatment response.

**Implication for the manuscript:** Present both compartments transparently. The basal spine data clearly support the sex-dependent narrative. The apical data, while showing a trend in the same direction, do not reach significance and should be discussed honestly.

## D.4 Converging evidence across endpoints

The consistent pattern across three independent endpoints strengthens the overall conclusion:

| Endpoint | Males effect size | Females effect size | Male p-value | Female p-value |
|----------|------------------|--------------------:|----------:|------------:|
| NOR (Cohen's d) | 1.26 (large) | 0.61 (medium) | 0.013 | 0.152 |
| Basal spines (Cliff's δ) | 0.49 (large) | 0.23 (small) | 0.008 | 0.216 |
| 3-NT (Cohen's d) | 2.75 (large) | 0.70 (medium) | 0.0002 | 0.223 |

The systematic attenuation in females — never a reversal, never an equivalent effect — argues against random variability and points to a shared underlying mechanism. Crucially, the cognitive endpoint (NORT), which is not confounded by baseline sex differences (both sexes show impairment relative to WT), shows the same pattern as the molecular and morphological markers, lending coherence to the interpretation.

This convergence sets the stage for investigating the mechanistic basis of sex-dependent PEA resistance. Candidate mechanisms include differential PEA metabolism (FAAH/NAAA activity), sex differences in PPAR-α expression or signaling, and estrogen-mediated modulation of neuroinflammatory pathways.
