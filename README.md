# Repository with the code to reproduce all the analyses of the paper _High resolution multi-scale profiling of embryonic germ cell-like cells derivation reveals pluripotent state transitions in humans_

**Sarah Stucchi, Lessly P. Sepulveda-Rincon, Camille Dion, Gaja Matassa, Alessia Valenti, Cristina Cheroni, Alessandro Vitriolo, Filippo Prazzoli, George Young, Marco Tullio Rigoli, Martina Ciprietti, Benedetta Muda, Zoe Heckhausen, Petra Hajkova, Nicolò Caporale, Giuseppe Testa, Harry G. Leitch**

<p align="center">
    <img src="https://www.biorxiv.org/content/biorxiv/early/2025/01/14/2025.01.14.632914/F1.large.jpg" width="550">
</p>


This repository contains all the code used to analyze the single-cell (scRNAseq) and DNA methylation (EMseq) data for the bioRxiv paper [High resolution multi-scale profiling of embryonic germ cell-like cells derivation reveals pluripotent state transitions in humans](https://doi.org/10.1101/2025.01.14.632914).

Docker images: 
- for most of the scRNAseq analyses (notebooks in [scRNASeq folder](scRNASeq) from 00 to 05) can be retrieved via `docker pull alessiavalenti/transgenerationalhub:transgenerational-1.1.7`.
- for the analyses of DNA methylation data (Rmarkdown in [EMSeq folder](EMSeq)) can be retrieved via `docker pull testalab/downstream:Transgenerational-1.1.5`.

An html version of the notebooks is accessible [here](https://GiuseppeTestaLab.github.io/EGCLC_paper_release/).
