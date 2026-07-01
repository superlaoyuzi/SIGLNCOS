# SIGLNCOS: Single-cell lncRNA Regulatory Landscape Analysis Framework

[![DOI](https://img.shields.io/badge/DOI-10.xxxx/xxxxxx-blue)](https://doi.org/10.xxxx/xxxxxx)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.2.1-276DC3)](https://www.r-project.org/)

SIGLNCOS is a computational framework designed to systematically investigate long non-coding RNA (lncRNA) regulation at single-cell resolution. It enables the identification of cellular signature lncRNAs (clncRNAs) and the characterization of their functional and co-regulatory relationships across diverse cancer types.

**Associated Paper**: *A Novel Framework for Identifying Cellular Signature LncRNAs and Their Functional Coregulatory Networks with Prognostic Implications at Single-Cell Resolution* (Journal of Advanced Research, 2026)

---

## Table of Contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Key Results](#key-results)
- [Installation](#installation)
- [Usage](#usage)
- [Data Availability](#data-availability)
- [Web Platform](#web-platform)
- [Repository Structure](#repository-structure)
- [Citation](#citation)
- [Contact](#contact)
- [License](#license)

---

## Overview

Long non-coding RNAs (lncRNAs) play important roles in tumor initiation, progression, and the tumor microenvironment. However, their systematic characterization across multiple cancer types at single-cell resolution remains limited.

To address this, SIGLNCOS provides a unified analytical workflow to construct pan-cancer single-cell lncRNA regulatory landscapes and identify functionally relevant lncRNA programs.

---

## Workflow

SIGLNCOS consists of four main analytical modules:

1. **Construction of lncRNA expression landscape**  
   Visualization and organization of lncRNA expression patterns across single cells from 31 cancer types (524,952 cells total).

2. **Identification of cellular signature lncRNAs (clncRNAs)**  
   Detection of lncRNAs with cell-type-specific expression patterns using differential expression analysis (log2FC > 0.25, adjusted p < 0.05, min.pct > 0.1).

3. **Functional characterization of clncRNAs**  
   Functional enrichment analysis (GO, KEGG, and immune-related gene sets) via GSEA to infer biological roles and pathway associations.

4. **Inference of co-regulatory lncRNA pairs**  
   Identification of lncRNA-lncRNA regulatory modules based on shared functional regulation patterns (Jaccard index > 0.5, |S| > 1).

Additionally, prognostic analyses (Cox regression and Kaplan-Meier) are performed to evaluate the clinical relevance of identified lncRNAs and co-regulatory networks.

---

## Key Results

Application of SIGLNCOS across 524,952 single cells from 31 cancer types led to the identification of:

- **667 cellular signature lncRNAs (clncRNAs)** associated with immune regulation, tumor progression, and cellular differentiation
- Functional evidence suggesting that **NEAT1** influences macrophage differentiation and cell fate decisions across cancer types (M1/M2 polarization analysis)
- **3,381 co-regulatory clncRNA pairs** involved in key cancer-related biological pathways

Co-regulatory analysis revealed complex lncRNA interaction networks that are particularly active during tumor evolution. For example, distinct co-lncRNA modules were identified in the transition from papillary thyroid carcinoma (PTC) to anaplastic thyroid carcinoma (ATC), suggesting coordinated regulatory mechanisms that cannot be explained by single lncRNAs alone.

Importantly, co-lncRNA network-based signatures showed improved performance in patient outcome prediction compared to individual lncRNAs (e.g., NR2F1-AS1 co-regulatory networks in BLCA fibroblasts, HR increased from 1.5 to 2.0), highlighting their potential clinical utility.

---

## Installation

### Option 1: Conda (Recommended)

```bash
# Clone the repository
git clone https://github.com/superlaoyuzi/SIGLNCOS.git
cd SIGLNCOS

# Create and activate conda environment
conda env create -f environment.yml
conda activate siglncos

# Install additional packages (if needed)
R -e "remotes::install_github('junjunlab/ClusterGVis')"
```

### Option 2: Docker

```bash
docker pull rocker/tidyverse:4.2.1
docker run -v $(pwd):/workspace -it rocker/tidyverse:4.2.1 R
```

### Option 3: Manual R Installation

```r
# Install core dependencies
install.packages(c('Seurat', 'harmony', 'igraph', 'ggplot2', 
                   'pheatmap', 'clusterProfiler', 'survival', 'survminer'))
BiocManager::install('org.Hs.eg.db')
```

---

## Usage

### Run the full pipeline

```bash
Rscript 01_core_scripts/01_SIGLNCOS_pipeline.R
```

### Run individual modules

```r
# Load functions
source("02_functions/02_RandomWalk.R")
source("02_functions/02_density_plot.R")

# Step 1: Identify clncRNAs
# Step 2: Functional enrichment
# Step 3: Co-regulatory network inference
# Step 4: Prognostic analysis
```

### Expected runtime
- Full analysis (31 cancer types): ~4-6 hours on HPC
- Single cancer type: ~30-60 minutes

---

## Data Availability

### Raw Data Sources
- scRNA-seq data: GEO/ArrayExpress (see Supplementary Table 1 for accession numbers)
- Bulk RNA-seq data: TCGA (GBM, HNSCC, BC, UCEC, CRC, OV)

### Processed Data
Download from: https://figshare.com/... (DOI: 10.xxxx/xxxxxx)

| File | Description |
|------|-------------|
| `all_colncRNA.txt` | 3,381 co-lncRNA pairs |
| `alldelncRNA3.txt` | 667 clncRNAs with annotations |
| `all_lnc_function.txt` | Functional annotations (GO/KEGG/Immune) |

### External Resources
| Resource | Version | Download Date |
|----------|---------|---------------|
| GENCODE | Release 43 (GRCh38.p13) | 2023 |
| ENCORI | v3.0 | 2023 |
| MSigDB | v2022.1 | 2023 |

---

## Web Platform

An interactive visualization platform, **sCo-Lnc**, is available for browsing, searching, and downloading:

**URL**: https://bio-bigdata.hrbmu.edu.cn/scolnc

### Features
- Browse by cancer type or lncRNA
- Search by gene symbol, Ensembl ID, or cell type
- Pan-cancer analysis views
- Download all results (node/edge tables)
- Online analysis tools (clustering, co-regulation, hallmark enrichment)

### Maintenance
- Regular updates every 1-2 years
- Minimum 5-year maintenance commitment
- Contact: scolnc@hrbmu.edu.cn for issues or feedback

---

## Repository Structure

```
SIGLNCOS/
├── 01_core_scripts/
│   └── 01_SIGLNCOS_pipeline.R    # Main analysis pipeline
├── 02_functions/
│   ├── 02_RandomWalk.R           # Random walk algorithm
│   ├── 02_density_plot.R         # Density plot functions
│   └── 02_odds_ratio.R           # Odds ratio calculation
├── 03_validation/                # Sensitivity analyses
├── 04_visualization/             # Optional visualization scripts
├── data/                         # Data directory (user-provided)
├── results/                      # Output directory
├── environment.yml               # Conda environment
├── LICENSE                       # MIT License
└── README.md                     # This file
```

---

## Citation

If you use SIGLNCOS or the sCo-Lnc database in your research, please cite:

> Zhang Y, et al. (2026). A Novel Framework for Identifying Cellular Signature LncRNAs and Their Functional Coregulatory Networks with Prognostic Implications at Single-Cell Resolution. *Journal of Advanced Research*. DOI: 10.xxxx/xxxxxx

---

## Contact

**Corresponding Author**: Prof. Yunpeng Zhang  
**Email**: yunpengzhang@hrbmu.edu.cn  
**Laboratory**: Bioinformatics Lab, Harbin Medical University  
**Database Support**: scolnc@hrbmu.edu.cn

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Last Updated**: 2026-07-01
```

---