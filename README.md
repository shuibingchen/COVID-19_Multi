# COVID-19_Multi

This repository contains the codes necessary to perform the analysis in our
manuscript ["A Multi-Organoid Platform Identifies CIART as a Key Factor for
SARS-CoV-2 Infection"](https://www.nature.com/articles/s41556-023-01095-y), as
described in the methods and main text.

### Input data

The single cell RNA-seq data were generated with the 10X Chromium and
pre-processed using 10X cellranger pipeline. The data are available in the
GEO database with accession#
[GSE202964](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202964).

The bulk RNA-seq data are available in the GEO database with accession#
[GSE202963](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202963).

The ATAC-seq data are available in the GEO database with accession#
[GSE202965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202965).

The CUT&RUN data are available in the GEO database with accession#
[GSE202966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202966).

### Requirements

The following R packages are needed

- Seurat
- scran
- scater
- batchelor
- ggplot2
- pheatmap
- RColorBrewer
- DiffBind
- dplyr
- tidyr
- tibble
- ChIPseeker
- DESeq2
- BiocParallel
- vsn
- plotly
- htmlwidgets
- karyoploteR

The following python packages are needed

- pandas

