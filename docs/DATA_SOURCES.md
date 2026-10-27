# Public inputs and snapshot verification

All file names below are relative to the project root. `config/inputs.csv` is the machine-readable contract; `config/input_sha256.tsv` records the manuscript input snapshots and the canonical CTD pair file prepared for this release. Input files are deliberately excluded from the code release. SHA-256 digests identify the exact input snapshots.

## TCGA expression and annotations

Obtain the TCGA-THCA STAR-count matrix from the [UCSC Xena GDC hub](https://xenabrowser.net/datapages/?dataset=TCGA-THCA.star_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net). Save the compressed TSV as `data/tcga/TCGA-THCA.star_counts.tsv.gz`. The expected input is gene by sample, with Ensembl identifiers in the first column and **log2(count + 1)** values. Raw integer counts, TPM or FPKM are not interchangeable with this input. The manuscript snapshot contains 572 tissue columns.

Download the complete HGNC tab-separated file from the [official HGNC download page](https://hgnc.genenames.org/download/) and save it as `data/reference/hgnc_complete_set.txt`. Required fields are `symbol` and `ensembl_gene_id`. For exact reproduction, retain the input snapshot matching the recorded fingerprint; the current download is mutable.

## Clinical data and mutations

Use the [cBioPortal THCA PanCancer Atlas study](https://www.cbioportal.org/study/summary?id=thca_tcga_pan_can_atlas_2018). Download its patient clinical file as `data/tcga/data_clinical_patient.txt`. Required fields are `PATIENT_ID`, `AGE`, `SEX`, `AJCC_PATHOLOGIC_TUMOR_STAGE`, `PATH_T_STAGE`, `PATH_N_STAGE` and `ICD_O_3_HISTOLOGY`. Public TCGA patient barcodes are used for matching. Selected public participant-level clinical and score fields are included in the curated plotting tables under `data/figure_tables`; raw clinical input files and private patient identities are not bundled. Only the listed curated public plotting inputs are bundled; obtain the raw clinical inputs separately.

Save the JSON string array from the [sequenced sample list](https://www.cbioportal.org/api/sample-lists/thca_tcga_pan_can_atlas_2018_sequenced/sample-ids) as `data/tcga/mutations/sequenced_samples.json`. Save the detailed mutation-record array returned by the [mutation endpoint](https://www.cbioportal.org/api/molecular-profiles/thca_tcga_pan_can_atlas_2018_mutations/mutations?sampleListId=thca_tcga_pan_can_atlas_2018_sequenced&projection=DETAILED) as `data/tcga/mutations/thca_mutations.json`. The analysis expects `keyword` and `sampleId` fields. A changed endpoint schema must be adapted explicitly rather than treating missing mutation records as wild type.

From the original [TCGA THCA publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4243044/), obtain the supplementary workbook named `NIHMS633017-supplement-15.xlsx`; save it under `data/reference/`. The `THCA-TP (496)` worksheet supplies BRS (`BRAF_RAF_score`), TDS (`differentiation_score`), `ERK_score` and `purity`. Only those fields and sample identifiers are read. The workbook is not redistributed.

## GEO

Download GPL570 platform annotation from [GEO GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) as `data/geo/GPL570.annot.gz`. Use the GEO annotation table with the `Gene symbol` field and its standard 27-line preamble. A raw manufacturer annotation with a different schema is not interchangeable.

Download compressed series matrix files from each accession's GEO page. Keep the original GEO header and expression scale. Save them using these exact file names:

| Accession | Relative file | Role |
|---|---|---|
| [GSE33630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33630) | data/geo/GSE33630_series_matrix.txt.gz | PTC and ATC versus normal |
| [GSE35570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35570) | data/geo/GSE35570_matrix.txt.gz | PTC after overlap removal |
| [GSE29265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29265) | data/geo/GSE29265_matrix.txt.gz | Supportive PTC assessment |
| [GSE60542](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60542) | data/geo/GSE60542_matrix.txt.gz | PTC after tissue selection and aggregation |
| [GSE65144](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65144) | data/geo/GSE65144_series_matrix.txt.gz | ATC transfer; MAS5 values log2 transformed before probe averaging |
| [GSE82208](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82208) | data/geo/GSE82208_series_matrix.txt.gz | FTC versus FTA transfer |

## CTD and online annotations

Obtain the chemical–gene interaction TSV from the [CTD download page](https://ctdbase.org/downloads/), preserving it as `data/ctd/CTD_chem_gene_ixns.tsv.gz`. Run `python tools/prepare_ctd.py` to create `data/ctd/ctd_human_chem_gene.tsv`: unique ChemicalName/GeneSymbol pairs for OrganismID 9606, sorted bytewise. The TSV contains two columns without a header. Its 829,489 pairs match the original analyzed pair membership; bytewise sorting replaces the original locale-dependent sorting. `config/ctd_null_gene_universe.txt` explicitly preserves the original 28,629-gene sampling order, so the fixed-seed null is reproducible across preparation locales. Step 15 stops if that gene universe changes. The chemical database is mutable; new snapshots require an explicit analysis update.

MSigDB sets are obtained through `msigdbr`, and GO annotation uses `org.Hs.eg.db`. KEGG pathway and STRING network inputs use the small frozen snapshots supplied under `data/annotation_snapshots`; the source manifest records their provenance. Other resource contents are not bundled. Dataset and annotation terms remain with their providers.

## Checking inputs

`python tools/check_inputs.py` distinguishes missing required files, absent optional files and changed snapshots using SHA-256. In `config/inputs.csv`, an empty `steps` field marks an optional or historical resource. The recorded doublet-call cache may be omitted for recomputation. Historical spatial geometry files retain their fingerprints but are not required by the current figure stages, which read the included frozen coordinate tables. A newer snapshot can be analyzed, but it should be identified as a new analysis rather than claimed to reproduce the manuscript exactly. The runner additionally records input MD5 checksums and package versions under the excluded `outputs/` directory.


Single-cell, spatial and composition inputs: raw UMI matrices from GSE193581, Visium matrices/features from GSE250521, and the upstream MCP-counter marker resource. Exact filenames and URLs are listed in config/inputs.csv. Place public data under the relative paths in that manifest. Step 18 generates broad cell-type annotations from the public count matrices.


## Reference snapshots

The included `data/annotation_snapshots/source_manifest.json` records public STRING v12.5 and KEGG URLs, verification on 19 October 2026 and SHA-256 digests. The GPL570 annotation header is dated 9 August 2016. The exact historical HGNC input is identified by digest in Supplementary Data 1; its original download date was not retained. MSigDB release 2026.1 was used with msigdbr 26.1.1 and org.Hs.eg.db 3.23.1. Replacing an input snapshot can change results.

The discovery set is the retained TCGA-THCA primary-tumour set. Clinical sample labels and ICD-O-3/GDC diagnoses are separate source fields; seven code-excluded records are documented in the Supplementary Data 1 histology provenance ledger. Eight primary tumours lack matched clinical/histology records. They remain in discovery and the seven-record sensitivity run.


## Epithelial-lineage analysis

Modules 32–33 review author-provided identities, apply the recorded per-sample doublet calls or recompute them, aggregate raw counts and fit independent-donor pseudobulk comparisons. Module 34 writes Figures 5–7 and Supplementary Figures S6–S7 after the prerequisite plotting modules. Spatial maps use retained spot coordinates on a neutral background. Historical histology images and spatial geometry files have empty `steps` fields in `config/inputs.csv`; they retain their fingerprints but are optional for the current workflow. Copy the `singlecell_lineage` directory from Supplementary Data 1 into `data/singlecell_lineage` in the analysis workspace. Obtain the GEO author annotation file `GSE193581_celltype_annotation.txt.gz` under `data/geo/GSE193581`. Module 18 supplies the processed sparse-count object. The supplied metadata identifies prior treatment and the two matched tissue pairs; no sample identity is inferred from outcome-associated expression.

For a fast reproduction of inference from the frozen pseudobulk sums, run:

```text
Rscript singlecell_lineage/reproduce_from_frozen.R path/to/Supplementary_Data_1/singlecell_lineage path/to/reproduced_lineage
```

This entry point checks all reported four-gene effect estimates, moderated intervals and P values against the supplied files. It reproduces inference from count sums; it does not repeat cell filtering or doublet classification. The primary tissue-labelled epithelial-lineage comparison (Q1) and exploratory malignant–normal follicular comparison (Q2) use three normal and four PTC specimens from seven independent donors. Alternative donor choices, original eight-specimen descriptive fits, donor omissions, the 100-cell threshold, treated material and the exploratory 30-cell identity unit are all retained. The complete eight-specimen fit is an unblocked descriptive comparator. Repeated-unit sensitivity models use two correlation rounds with block-aware voom between them. Confidence intervals and P values share the empirical Bayes model.

To force fresh doublet calling, omit the input `data/singlecell_lineage/dbl_calls_cache.csv.gz` and the generated cache under `outputs/analysis/extended/singlecell_lineage/01_annotation`. An existing generated cache is reused; when it is absent, a supplied input cache is copied before recomputation is considered.
