# A deep cellular atlas of the human ventral substantia nigra in Parkinson’s identifies a genetic and molecular overlap with insulin resistance.

## Summary

This repository contains the code to reproduce figures and analysis presented in the manuscript: A deep cellular atlas of the human ventral substantia nigra in Parkinson’s identifies a genetic and molecular overlap with insulin resistance, Volpato et al. (2026)

## Data availability

Raw and processed gene expression data has been deposited in GEO under the accession number:

## Repository Structure

- The scripts/ directory includes all analyses discussed in the manuscript

```
deep_snRNAseqAtlas_humanSN/
├── README.md
├── LICENSE
├── scripts/
│   ├── 01_data_preprocessing.R
│   ├── 02_cell_type_annotation.R
│   ├── 03_subtype_analysis.R
│   ├── 04_DEG_analysis.R
│   ├── 05_DTU_analysis.R
│   ├── 06_pathway_PPI_analysis.R
│   ├── 07_pseudotime_analysis.R
│   ├── 08_cell_proportions.R
│   ├── 09_cellcommunication_analysis.R
│   ├── 10_coexpression_network.R
│   ├── 11_SCENIC.R
│   ├── cellchat_utils.R
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
│   ├── Figure_3.Rmd
│   ├── Figure_4.Rmd
│   ├── Figure_5.Rmd
│   └── Figure_6.Rmd

```

## Software
| **Software**                        | **Version(s)** | **Notes**                                                          |
| ----------------------------------- | -------------- | -------------------------------------------------------------------------------------------------------------------------------- | 
| MAPPA pipeline (Takara)                             | v 1.0   | Used for raw sequencing data processing.                                       |
| Salmon                            | v.1.4.0        |  Used for transcript quantification.                                |
| MAGMA                            | v.1.08        |  Used for gene-set association analyses.                                |
| gcta                            | v1.92.2beta        |  Used for GWAS conditional analysis.                                |
| ggplot2                             | v.3.5.0        |  Used for data visualization in R.                                  |
| propeller                           | v.1.0.0        |  Used for differential cell type proportion analysis in R.               |
| Seurat                              | v.4.1.3        |  Used for single-nucleus RNA-seq analysis in R.          |
| topGO                             | v.2.61.1        |  Used for Gene Ontology enrichment analysis in R.                                  |
| CellChat                             | v.2        |  Used for cell-cell communication analysis in R.                                  |
| rrvgo                             | v.3.2        |  Used for reduction of Gene Ontology enriched terms in R.                                  |
| igraph                             | v.2.1.4        |  Used for gene network module clustering in R.                                  |
| scProportionTest                             | [link](https://github.com/rpolicastro/scProportionTest)        |  Used for differential cell type proportion analysis in R.                                  |
| slingshot                             | v.1.8        |  Used for pseudotemporal ordering analysis in R.                                  |
| switchde                             | v.1.35.0        |  Used for identification of differentially expressed genes along pseudotime trajectory in R. |
| tradeSeq                             | v.1.5.10        |  Used for identification of differentially expressed genes along pseudotime trajectory in R  |
| BigScale                             | v.2.0        |  Used for gene network analysis in R.                                  |
| SCENIC                             | v.1.3.1        |  Used for TF regulons analysis in R.                                  |
| Fishpond                             | v.2.16        |  Used for transcript analyses in R.                                  |
| R Project for Statistical Computing | v.4.3.0        |  Used for statistical computing and data analysis.                  |
| PLINK                               | v.2.0   |  Used for genetic and association analyses as described here https://github.com/csandorfr/pd-t2d-meta-gwas                         |
| METAL                               | build 2011-03-25   |  Used for fixed-effect inverse-variance meta-analysis as described here https://github.com/csandorfr/pd-t2d-meta-gwas                         |

