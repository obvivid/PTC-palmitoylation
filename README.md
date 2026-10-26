# Palmitoylation-related expression signature in papillary thyroid carcinoma

Analysis code accompanying *A palmitoylation-related gene expression signature distinguishes papillary thyroid carcinoma across four GEO datasets*. It covers TCGA discovery, patient-grouped nested validation, GEO evaluation, molecular and clinical associations, single-cell and spatial analysis, exploratory individual-profile scoring, and the manuscript figures.

## Requirements and data

The manuscript used R 4.6.0 and the package versions recorded in `config/package_versions.tsv`. Python 3.10 or newer, using the standard library, is needed for input preparation and checking. `install_dependencies.R` installs missing available packages; it does not pin them to the recorded versions.

Obtain the public inputs following [DATA_SOURCES.md](docs/DATA_SOURCES.md) and `config/inputs.csv`, and save them under the listed relative `data/` paths. For the epithelial-lineage analysis, also copy `singlecell_lineage` from Supplementary Data 1 to `data/singlecell_lineage`. Required raw assay matrices are obtained separately. The included `data/figure_tables` and `data/annotation_snapshots` contain the small plotting inputs and reference snapshots used by the scripts. Public-study specimen identifiers and selected public clinical fields are retained where needed for traceability.

## Run the analyses

Run from this repository folder:

```text
Rscript install_dependencies.R
python tools/prepare_ctd.py
python tools/check_inputs.py
Rscript run_analysis.R --dry-run
Rscript run_analysis.R
```

Prepare the CTD input before `prepare_ctd.py`, as described in `docs/DATA_SOURCES.md`. The dry run lists the analysis steps and missing input files without performing analyses. Results are written to `outputs/`.

The 34 modules and their prerequisites are listed in `config/steps.csv`. To run a subset, use `Rscript run_analysis.R --steps=01,02,03`; prerequisite outputs must already exist. Nested model comparisons and single-cell preprocessing take longer than the other modules. Analytical definitions are in [METHODS.md](docs/METHODS.md), application details in [APPLICATION.md](docs/APPLICATION.md), and the final figure mapping in [PUBLICATION_FIGURES.md](docs/PUBLICATION_FIGURES.md). Modules 30 and 34 write the coordinated and final revised figures, respectively.

## Scoring examples

The equal-weight expression contrast and exploratory learned rule use separate frozen development parameters:

```text
Rscript single_sample/score_cli.R single_sample/example_four_gene_input.csv result_contrast.csv contrast youden
Rscript application_models/score_cli.R application_models/example_four_gene_input.csv result_weighted.csv
```

Inputs are deposited, probe-averaged GPL570 gene-level log2 expression. These examples demonstrate execution. The TCGA signature, target-cohort-standardised GEO score and fixed individual-profile rules are distinct definitions; their scope is described in the manuscript and application documentation.

To reproduce donor-level inference from the supplied frozen pseudobulk counts:

```text
Rscript singlecell_lineage/reproduce_from_frozen.R path/to/Supplementary_Data_1/singlecell_lineage path/to/reproduced_lineage
```

This entry point reproduces inference from count sums. Modules 18 and 32 perform the upstream cell processing. Doublet-call recomputation requires omitting the supplied input cache and the corresponding generated cache; see `docs/DATA_SOURCES.md`.

## Analysis checks

Optional checks of patient-level fold assignment, clinical-field parsing and CTD species filtering use small synthetic inputs:

```text
Rscript tests/test_portability.R
python -m unittest discover -s tests
```

Input snapshots are checked by `tools/check_inputs.py`. Result-level checks are retained in the analytical modules, including the frozen scoring and mixing checks under `application_models/`.

## Licence and citation

The original analysis code is distributed under the [MIT License](LICENSE). Public datasets, annotation snapshots and software dependencies retain their providers' terms. Citation metadata are in `CITATION.cff`; please cite the associated paper when using this work.
