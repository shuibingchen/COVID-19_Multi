# step_1.ambient_RNA_removal.R
# step 1. removal of ambient RNA (background noise) in UMI counts matrix with SoupX
# 

library(SoupX)
library(DropletUtils)
library(tidyverse)
library(Matrix)

workdir <- "."
srcdir <- file.path(workdir, "src", "scRNAseq")
umidir <- file.path(workdir, "data", "UMI_counts")

# load functions
source(file.path(srcdir, "my_functions.R"))

# set a random seed
set.seed(98)

# perform ambient RNA removal for each sample
for (sid in c("CIART_Mock","CIART_SARS_CoV_2","WT_Mock","WT_SARS_CoV_2")){
  vtable <- my.soupx(umidir, sid)
}
