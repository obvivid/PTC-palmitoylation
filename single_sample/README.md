# Fixed single-sample scoring

Run the scoring entry point with a CSV containing gene and log2_expression columns:

```text
Rscript single_sample/score_cli.R single_sample/example_four_gene_input.csv result.csv contrast youden
```

The expression contrast requires the four genes only. The rank option requires the complete frozen gene background. Inputs use probe-averaged GPL570 gene-level log2 expression as specified in Supplementary Methods SM11. The primary contrast threshold is -0.363364815791667, and scores at or above it are positive. Parameters are fixed in the supplied text files; scoring uses no target-cohort statistics. The example is a four-gene excerpt of a public development profile and demonstrates execution only.

This module implements an exploratory expression scoring procedure. Reproduction modules 25 and 26 regenerate its retrospective evaluation and Supplementary Figure S9. Synthetic mixing weights describe expression signals. Raw-sample processing and clinical interpretation follow the scope reported in the manuscript.
