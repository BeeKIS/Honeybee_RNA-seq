## Honeybee RNA-seq Analysis

### Introduction
This section provides an overview of the data and code files used for RNA-seq analysis of honeybees.

### Data:
- `counts_data.csv`: Count matrix used for Differential Gene Expression (DGE) analysis with DESeq2
- `sample_info.csv`: Metadata on hygienic behavior between samples
- `allOE_genes.txt`: List of all overexpressed genes used for functional analysis
- `sigOE_genes.txt`: List of significantly overexpressed genes used for functional analysis

### `DESeq2_analysis.R`:
R script for conducting Differential Gene Expression (DGE) analysis using DESeq2.

### `Gene_annotation.R`:
R script for gene annotation.

### `enrichGO_analysis.R`:
R script for functional analysis.
