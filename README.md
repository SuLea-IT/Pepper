# Pepper Code Repository

This repository stores analysis scripts for our pepper (Capsicum) single-nucleus and spatial transcriptomics study.

## Publication

**Title:** *Laminar patterning transcription factors orchestrate spatial metabolite partitioning in Capsicum fruit*

## Abstract

Chili pepper (*Capsicum annuum* L.) produces specialized metabolites, notably the pungent capsaicin and the red capsanthin. Although their biosynthetic pathways are well characterized, the cellular architecture that underpins spatial regulation remains unclear. Here, we present a spatiotemporal single-nucleus atlas of pepper development, integrating single-nucleus RNA sequencing and spatial transcriptomics, profiling 332,468 high-quality cells from 57 samples spanning seedlings to mature fruits. This resource reveals a multilayered organization and precisely maps metabolic genes to defined cell types and spatial regions. We further identify laminar patterning transcription factors (LPTFs), including WRKY6, ZAT10, and BTF3, whose layer-specific expression correlates with localized capsanthin accumulation. Our work establishes a framework for dissecting laminar control of specialized metabolism and provides a valuable reference for comparative studies across species. The atlas is openly accessible at http://Pepper-Cell-Atlas.com.

## Required R packages

The following packages are required:

`Seurat`, `ggplot2`, `patchwork`, `dplyr`, `here`, `tidyverse`, `viridis`, `lattice`, `reshape2`, `cowplot`, `Matrix`, `Matrix.utils`, `edgeR`, `S4Vectors`, `SingleCellExperiment`, `pheatmap`, `apeglm`, `png`, `DESeq2`, `RColorBrewer`, `data.table`

## Repository structure

- `Script/`: analysis scripts for clustering, integration, projection, marker analysis, pseudotime, and visualization.
