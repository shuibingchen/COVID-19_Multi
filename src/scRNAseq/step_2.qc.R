# step_2.qc.R
# step 2: perform QC based on the corrected UMI counts table
# 

library(Seurat)
library(magrittr)
library(future)

workdir <- "."
srcdir <- file.path(workdir, "src", "scRNAseq")
sourcedir <- file.path(workdir, "data")
umidir <- file.path(sourcedir, "UMI_counts")
figdir <- file.path(workdir, "result", "scRNAseq", "figure")
infodir <- file.path(workdir, "result", "scRNAseq", "info")

# project name
project <- "shuibing"

# pattern for defining mitochondrial/ribosomal/virus genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"
virus.pattern <- "^CoV2"

# set a random seed
set.seed(98)

# load functions
source(file.path(srcdir, "my_functions.R"))

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# sample info
sample.info <- data.frame(SeqName=c("CIART_Mock","CIART_SARS_CoV_2","WT_Mock","WT_SARS_CoV_2"), 
                          Name=c("CIART_Mock","CIART_SARS-CoV-2","WT_Mock","WT_SARS-CoV-2"), 
                          Type=c('CIART-/-','CIART-/-','WT','WT'), 
                          Virus=c('Mock','SARS-CoV-2','Mock','SARS-CoV-2'))
rownames(sample.info) <- c("XM1","XM2","XM3","XM4")

# load raw UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  pid <- rownames(sample.info)[k]
  sid <- sample.info$SeqName[k]
  raw.counts.list[[k]] <- my.Read10X(file.path(umidir, sid, "soupx", "matrix"), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix(raw.counts.list)

# Initialize the Seurat object with the raw (non-normalized data).
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
                                   names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)
panc.initial[["percent.virus"]] <- PercentageFeatureSet(panc.initial, pattern=virus.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc.initial@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.initial@meta.data[,"orig.ident"])])
}
panc.initial %<>% AddMetaData(metadata=tmeta)

# calculate nGenes/nUMIs after excluding virus genes
non.virus.genes <- grep(virus.pattern, rownames(raw.counts.all), value=T, invert=T)
# Initialize a temporary Seurat object excluding virus genes.
temp.seurat <- CreateSeuratObject(counts=raw.counts.all[non.virus.genes, ], project=project, assay="RNA", min.cells=0, min.features=0, 
                                  names.field=1, names.delim="_", meta.data=NULL)
nonvirus.stats <- FetchData(temp.seurat, vars=c('nCount_RNA','nFeature_RNA')) %>% rename('nCount_RNA_nonVirus'='nCount_RNA','nFeature_RNA_nonVirus'='nFeature_RNA')
panc.initial %<>% AddMetaData(metadata=nonvirus.stats)

# perform cell filtering
# nGene > 500, nGene <= 8000, nUMI > 1000, nUMI <= 70000, percent.mito < 15%
panc.initial %<>% subset(subset=nFeature_RNA_nonVirus > 500 & nFeature_RNA_nonVirus <= 8000 & nCount_RNA_nonVirus > 1000 & nCount_RNA_nonVirus <= 70000 & percent.mito < 15)

# free space
rm(temp.seurat)
rm(raw.counts.list)
rm(raw.counts.all)
rm(tmeta)
rm(tdic)
gc()

# save Seurat object (initial)
saveRDS(panc.initial, file=file.path(infodir, "panc.initial.rds"))
