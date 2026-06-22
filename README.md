# SIGLNCOS: Single-cell lncRNA Regulatory Landscape Analysis Framework

SIGLNCOS is a computational framework designed to systematically investigate long non-coding RNA (lncRNA) regulation at single-cell resolution. It enables the identification of cellular signature lncRNAs and the characterization of their functional and co-regulatory relationships across diverse cancer types.

## Overview

Long non-coding RNAs (lncRNAs) play important roles in tumor initiation, progression, and the tumor microenvironment. However, their systematic characterization across multiple cancer types at single-cell resolution remains limited.

To address this, SIGLNCOS provides a unified analytical workflow to construct pan-cancer single-cell lncRNA regulatory landscapes and identify functionally relevant lncRNA programs.

## Workflow

SIGLNCOS consists of four main analytical modules:

1. **Construction of lncRNA expression landscape**  
   Visualization and organization of lncRNA expression patterns across single cells.

2. **Identification of cellular signature lncRNAs**  
   Detection of lncRNAs associated with specific cellular identities or states (clncRNAs).

3. **Functional characterization of clncRNAs**  
   Functional enrichment analysis to infer biological roles and pathway associations.

4. **Inference of co-regulatory lncRNA pairs**  
   Identification of lncRNA–lncRNA regulatory modules based on shared functional regulation patterns.

Additionally, prognostic analyses are performed to evaluate the clinical relevance of identified lncRNAs and co-regulatory networks.

## Key Results

Application of SIGLNCOS across 524,952 single cells from 31 cancer types led to the identification of:

- **667 cellular signature lncRNAs (clncRNAs)** associated with immune regulation, tumor progression, and cellular differentiation
- Functional evidence suggesting that **NEAT1** may influence macrophage differentiation and cell fate decisions
- **3,381 co-regulatory clncRNA pairs** involved in key cancer-related biological pathways

Co-regulatory analysis revealed complex lncRNA interaction networks that are particularly active during tumor evolution. For example, distinct co-lncRNA modules were identified in the transition from papillary thyroid carcinoma to anaplastic thyroid carcinoma, suggesting coordinated regulatory mechanisms that cannot be explained by single lncRNAs alone.

Importantly, co-lncRNA network-based signatures showed improved performance in patient outcome prediction compared to individual lncRNAs, highlighting their potential clinical utility.

## Conclusion

SIGLNCOS provides a systematic framework for constructing pan-cancer single-cell lncRNA regulatory landscapes. It enables the discovery of both individual and cooperative lncRNA regulatory mechanisms underlying tumor biology.

This work offers new insights into cancer-associated lncRNA functions and supports the development of potential biomarkers and therapeutic targets.

## Web Platform

An interactive visualization platform is available:

http://www.bio-bigdata.hrbmu.edu.cn/scolnc/

## Notes

This repository contains scripts and workflows used in the analysis described above, supporting reproducibility and further exploration of single-cell lncRNA regulation in cancer.
