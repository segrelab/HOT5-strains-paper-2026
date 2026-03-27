# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research repository for the Osborne et al. 2026 manuscript on whole genome sequencing and trait-based/metabolic analyses of five HOT5 bacterial strains: **B08, C03, B10, E06, and F03**. The repo is organized around figures and analyses, each in its own folder.

## Repository Structure

- **`genome-assembly/`** — Scripts (by E.F. and C.H.) for assembly with Flye and polishing with Racon.
- **`Figure1-phylogenetic-tree/`** — R scripts and data files for the phylogenetic tree figure. Tree was generated in KBase using "Insert Genome Into SpeciesTree v2.2.0"; output Newick files (`new-tree.newick`, `new-tree-labels.newick`) are the inputs to the R scripts.

## Running R Scripts

No build system. Scripts are run directly in R or RStudio:

```r
# From within the Figure1-phylogenetic-tree/ directory:
source("Figure1.R")   # Main phylogenetic tree figure
source("eddies-71.R") # 71-species phylogenomic tree
source("SCGs-tree.R") # Single-copy gene tree (Rhodobacter)
```

**Required R package:** `ape` (phylogenetics)

**Important:** `Figure1.R` has a hardcoded `setwd("~/npsegre/melisa2026")`. Update this path to the actual `Figure1-phylogenetic-tree/` directory before running.

## Key Data Files (Figure1-phylogenetic-tree/)

| File | Description |
|------|-------------|
| `new-tree-labels.newick` | Primary input tree with KBase internal labels (used by Figure1.R) |
| `new-tree.newick` | Alternate Newick output from KBase SpeciesTree |
| `Bac_71_concatenated-proteins.fa` | Concatenated protein alignments for 71 bacterial species |
| `phylo-profile.db` | SQLite database of phylogenomic profile data |

## Analysis Workflow

```
Raw nanopore reads
    → Genome assembly (Flye) + polishing (Racon)
    → Genome annotation (KBase RAST/RASTtk)
    → Import into KBase narrative (narrative/206763)
    → Phylogenetic tree (KBase SpeciesTree v2.2.0)
    → Export Newick → R (ape) → Figure 1
```

The 5 HOT5 genomes plus reference genomes from Valiya Kaladi et al. 2026 (imported via NCBI accession numbers) are included in the tree.

## Planned Analyses (not yet in repo per README)

- Metabolic percolation analysis (D.B.)
- Trait-based analysis (from Zoccarato et al. 2022 methods)
- Raw data hosted on Dryad
