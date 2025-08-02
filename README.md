# In silico analysis of plasticity in the pancreatic lineages

This repository contains the code and data processing pipeline used in my master's dissertation for reconstructing cell lineages and analysing single-cell RNA-seq data from human pancreas development using CellTag, Seurat, and GEMLI.  
Raw data, metadata, and outputs are not available due to unpublished information.

## Analyses were performed using:

- [Seurat](https://github.com/satijalab/seurat) – scRNA-seq quality control, processing, and clustering  
- [CellTag](https://github.com/morris-lab/BiddyetalWorkflow.git) – CellTag-mediated lineage tracing  
- [GEMLI](https://github.com/UPSUTER/GEMLI) – in silico lineage inference  
- R version 4.4.1 on macOS Sonoma 14.2.1  

The complete list of R packages and versions can be found in the Methods section of the dissertation.

---

## Overview

This project investigates lineage plasticity during human pancreatic development using single-cell RNA-seq data. We benchmarked experimental and computational lineage tracing approaches, focusing on their complementarity in capturing developmental trajectories from induced pluripotent stem cells and fetal tissues.


## Computational Methods

### 1. CellTag Lineage Tracing

We implemented and modified the original [CellTag workflow](https://github.com/morris-lab/BiddyetalWorkflow) to accommodate ScaleBiosciences-style barcodes and quantified CellTag abundance using:

- Custom motif matching on BAM and FASTQ files  
- CellTag UMI parsing and whitelist filtering  
- Cell-by-cell barcode matrix construction  
- Clone inference based on barcode sharing  

### 2. Single-cell RNA-seq Analysis

Performed using the Seurat package:

- Doublet detection with `scDblFinder`  
- SCTransform normalisation and regression  
- Dimensionality reduction: PCA, UMAP, and t-SNE  
- Clustering and marker-based cell type annotation  
- Trajectory inference via diffusion maps and `Slingshot`  

### 3. GEMLI-based In Silico Lineage Reconstruction

Utilised [GEMLI](https://github.com/UPSUTER/GEMLI) to infer lineages based on gene expression memory:

- Performed on a clinical fetal pancreas dataset  
- Custom Singularity container created for HPC deployment  
- Lineages annotated and visualised based on cell types  
- Both symmetrical and asymmetrical lineages were identified with high confidence  

## Repository Contents

- `readmes/`: Scripts and markdown files showing the analysis workflows and their purposes.
  - `CellTag_BAM_processing/`: Scripts to extract and quantify CellTag motifs from aligned BAM files (e.g., regex filtering, CellTag parsing).
  - `CellTag_lineages/`: R scripts to assign CellTags to cells, generate CellTag-by-cell matrices, and infer clonal relationships based on barcode overlap.
  - `CellTag_seurats/`: Code to process and analyse Seurat objects from CellTagged scRNA-seq data, including QC, normalisation, clustering, and pseudotime analysis.
  - `Ma_seurats/`: Seurat analysis scripts for the public fetal pancreas dataset from Ma et al., used for in silico lineage tracing with GEMLI. 

## Citation

If you use any component of this pipeline or analysis, please cite:

> Ma, S. (2025). *In silico analysis of plasticity in the pancreatic lineages*. Master’s thesis, King’s College London.
