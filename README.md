# A deep cellular atlas of the human ventral substantia nigra in Parkinson’s identifies a genetic and molecular overlap with insulin resistance.

This repository contains the code to reproduce figures and analysis presented in the manuscript: A deep cellular atlas of the human ventral substantia nigra in Parkinson’s identifies a genetic and molecular overlap with insulin resistance, Volpato et al. (2026)

## Data availability:

Raw and processed gene expression data has been deposited in GEO under the accession number:

## Repository Orientation:

### The analysis/ directory includes all analyses discussed in the manuscript

```
deep_snRNAseqAtlas_humanSN/
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_data_preprocessing.R
│   ├── 02_cell_type_annotation.R
│   ├── 03_subtype_analysis.R
│   ├── 04_DEG_analysis.R
│   ├── 05_pathway_PPI_analysis.R
│   ├── 06_pseudotime_analysis.R
│   ├── 07_cell_proportions.R
│   └── utils.R                # Common helper functions
├── figures/
│   ├── Figure_1/
│   ├── Figure_2/
│   ├── Figure_3/
│   ├── Figure_4/
│   ├── Figure_5/
│   └── Figure_6/
├── notebooks/
│   ├── Figure_1.Rmd
│   ├── Figure_2.Rmd
│   └── Figure_3.Rmd
├── results/
│   ├── DEG_tables/
│   ├── MAGMA/
│   └── GO_PPI/
```
