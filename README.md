# In silico analysis of plasticity in the pancreatic lineages

This repository contains the code and data processing pipeline used in my master's dissertation for reconstructing cell lineages and analyzing single-cell RNA-seq data from human pancreas development using CellTag, Seurat, and GEMLI.

---

## Project Structure
scripts/ # scripts for each step of analysis
metadata/ # metadata summarising the analysis result

## Analyses were performed using:

- [Seurat](v5.3.0) – scRNA-seq quality control, processing, and clustering
- [CellTag] (https://github.com/morris-lab/BiddyetalWorkflow.git) - CellTag-mediated lineage tracing
- [GEMLI](v0.1.0) – in silico lineage inference  
- R version 4.4.1 on macOS Sonoma 14.2.1

The full list of R packages and versions can be found in Method section of the dissertation.
