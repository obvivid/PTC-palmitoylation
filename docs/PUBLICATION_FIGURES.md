# Figure reproduction guide

Step 30 prepares the coordinated plotting objects and complete cellular summaries. Step 34 writes the final Figures 5–7 and Supplementary Figures S6–S7. Curated public plotting tables supply the retained coordinates and annotations; figure rendering does not repeat feature selection or model training.

| Figure | Content |
| --- | --- |
| 1 | Candidate expression heat map and differential-expression volcano plot |
| 2 | Three-method overlap, fitted coefficients, expression, score, repeated nested AUCs and selection frequencies |
| 3 | External ROC curves, AUC confidence intervals and gene expression distributions |
| 4 | GO biological processes and STRING associations, with the four selected genes highlighted |
| 5 | Mutation status, published molecular axes, tumour purity and mutation-stratified tissue differences |
| 6 | Composition-adjusted coefficients, broad cell types, transcript detection and independent-donor epithelial-lineage expression |
| 7 | Spatial detection maps in the first numbered normal and PTC sections, and detection summaries for all 16 sections |
| S1 | Paired tumour–normal expression differences |
| S2 | KEGG enrichment |
| S3 | Chemical–gene association counts and random backgrounds |
| S4 | Clinicopathological score distributions |
| S5 | Sample-level thyrocyte detection and spatial section summaries |
| S6 | Donor-level expression, identity contrast and donor-omission sensitivity |
| S7a–S7d | Complete spatial detection maps for all four tissue groups |
| S8 | ATC comparisons and follicular carcinoma–adenoma discrimination |
| S9 | Fixed single-sample scoring and expression-signal mixing |
| S10 | Exploratory model comparison and weighted mixing |
| S11 | Canonical-marker display across all 22 clusters |

Figure 7a contains two rows (normal N-1 and PTC PTC-1) and two columns (ABHD17C and ZDHHC20). Main and supplementary spatial maps use a neutral background, with grey points for undetected transcripts and orange points for raw counts greater than zero. Maps within each section use identical coordinates, scale and point sizes. Counts and percentages describe the displayed section; section labels do not assign malignant-cell identity to individual spots. Figure 7b uses one point per section and black median bars.

Spatial figures use the supplied coordinates and detection states for 52,906 retained spots across 16 sections. Rendering Figure 7 and Supplementary Figure S7 requires no histology-image file. The geometry and detection summaries remain traceable to the public GEO study. No deposited tissue image is reproduced in these final spatial figures.

Gene and RNA symbols are italicised in expression plots. Protein labels in the STRING association network, family names and molecular-score names retain roman type. Mixed labels use ggtext to preserve gene italics alongside sample counts and ordinary text. The final figure files are produced by modules 30 and 31, with Figures 5–7 and Supplementary Figures S6–S7 subsequently written by module 34. This typography update does not alter analysis inputs, model parameters or numerical results.

The manuscript contains seven main figures and Table 1. PNGs are rendered at 300 dpi, TIFFs at 600 dpi and PDFs retain vector marks where supported. The supplied plotting tables contain public study identifiers, coordinates and summary values.
