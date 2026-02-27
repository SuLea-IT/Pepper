# Pepper Code Repository

This repository stores analysis scripts for our pepper (Capsicum) single-nucleus and spatial transcriptomics study.

## Publication

**Title:** *Laminar patterning transcription factors orchestrate spatial metabolite partitioning in Capsicum fruit*

## Abstract

Chili pepper (*Capsicum annuum* L.) produces specialized metabolites, notably the pungent capsaicin and the red capsanthin. Although their biosynthetic pathways are well characterized, the cellular architecture that underpins spatial regulation remains unclear. Here, we present a spatiotemporal single-nucleus atlas of pepper development, integrating single-nucleus RNA sequencing and spatial transcriptomics, profiling 332,468 high-quality cells from 57 samples spanning seedlings to mature fruits. This resource reveals a multilayered organization and precisely maps metabolic genes to defined cell types and spatial regions. We further identify laminar patterning transcription factors (LPTFs), including WRKY6, ZAT10, and BTF3, whose layer-specific expression correlates with localized capsanthin accumulation. Our work establishes a framework for dissecting laminar control of specialized metabolism and provides a valuable reference for comparative studies across species. The atlas is openly accessible at http://Pepper-Cell-Atlas.com.

## Required packages

The following are the packages currently used by scripts in this repository.

### R packages (used in `Script/*.R`)

`clustree`, `dplyr`, `future`, `ggplot2`, `ggpubr`, `harmony`, `Matrix`, `monocle`, `openxlsx`, `pheatmap`, `png`, `RColorBrewer`, `readr`, `readxl`, `reshape2`, `scales`, `sctransform`, `Seurat`, `SeuratDisk`, `stringr`, `tibble`, `tidyverse`, `viridis`

### Python package (used in `Script/merge_excel_files.py`)

`pandas`

## Repository structure

- `Script/`: analysis scripts for clustering, integration, projection, marker analysis, pseudotime, and visualization.
