# Analysis definitions

The detailed settings are in Supplementary Methods SM1–SM13 and the numbered modules.

| Analysis | Definition and implementation |
|---|---|
| Discovery | 505 TCGA-THCA primary tumours and 59 adjacent normals; sample identifiers matched to matrix columns; 31 predefined enzyme candidates, 29 expressed |
| Differential expression | Patient-blocked limma, consensus residual correlation, robust/trend empirical Bayes; BH across 17,697 rows |
| Selection | LASSO, top-ten random forest and linear SVM-RFE intersection; complete-data logistic fit; nominal patient-clustered HC1 coefficient intervals |
| Internal validation | Five repeats of ten patient-grouped outer folds; training-fold blocked DE, selection and fit; grouped inner LASSO; fixed mapping/universe |
| External score | ΣβZ after standardisation within each target comparison; no refit; patient bootstrap primary interval, DeLong secondary |
| External gene tests | OLS group contrast with patient-clustered HC1 covariance; t reference by patient count; BH across four genes per dataset |
| Clinical | Univariable rank or correlation tests; BH over eight variables; binary median-difference and age-correlation bootstrap intervals |
| Composition and mutation strata | OLS point estimates, patient-clustered HC1 covariance with cluster correction; reduced proxy checks; three contrasts per outcome BH-adjusted |
| Cellular | Maximum canonical-marker mean per cluster, verified in step 31; descriptive detection, positive-cell-count sample summaries and explicit zero-thyrocyte missingness |
| Spatial | Raw-count detection by section and pooled tissue group; retained deposited coordinates and image scales; section units |
| Annotation | GO/KEGG annotation context; MSigDB within-collection BH; STRING functional associations, score 400, no added genes |
| CTD | 1,000 random 25-gene sets; upper-tail empirical statistics; three-statistic BH |
| Survival | Four original source endpoints, median score groups, univariable Cox and log-rank; endpoint-specific missingness and four-endpoint BH |
| Application | Separate frozen C/R/L scoring definitions, training-only transformations/thresholds; retrospective transfer and descriptive matched expression-signal mixing |

`stable_wilcox_test` calculates the smaller direct one-sided tail twice for two-sided asymptotic tests. It avoids subtracting a rounded-to-one CDF. P values from this numerical fix and patient-clustered mean-model inference refer to distinct tests.

See `docs/PUBLICATION_FIGURES.md` for the final figure mapping and `docs/APPLICATION.md` for application modelling settings.


## Epithelial-lineage analysis

Modules 32–33 review author-provided identities, apply the recorded per-sample doublet calls or recompute them, aggregate raw counts and fit independent-donor pseudobulk comparisons. Module 34 writes the final revised Figures 5–7 and Supplementary Figures S6–S7 after the earlier plotting modules. Spatial maps use retained spot coordinates on a neutral background. Copy the `singlecell_lineage` directory from Supplementary Data 1 into `data/singlecell_lineage` in the analysis workspace. Obtain the GEO author annotation file `GSE193581_celltype_annotation.txt.gz` under `data/geo/GSE193581`. Module 18 supplies the processed sparse-count object. The supplied metadata identifies prior treatment and the two matched tissue pairs; no sample identity is inferred from outcome-associated expression.

For a fast reproduction of inference from the frozen pseudobulk sums, run:

```text
Rscript singlecell_lineage/reproduce_from_frozen.R path/to/Supplementary_Data_1/singlecell_lineage path/to/reproduced_lineage
```

This entry point checks all reported four-gene effect estimates, moderated intervals and P values against the supplied files. It reproduces inference from count sums; it does not repeat cell filtering or doublet classification. Main Q1 and exploratory Q2 use three normal and four PTC specimens from seven independent donors. Alternative donor choices, original eight-specimen descriptive fits, donor omissions, the 100-cell threshold, treated material and the exploratory 30-cell identity unit are all retained. The complete eight-specimen fit is an unblocked descriptive comparator. Repeated-unit sensitivity models use two correlation rounds with block-aware voom between them. Confidence intervals and P values share the empirical Bayes model.
