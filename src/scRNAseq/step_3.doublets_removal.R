# step_3.doublets_removal.R
# step 3: detect doublets with DoubletFinder
# 

library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(Matrix)
library(magrittr)

workdir <- "."
srcdir <- file.path(workdir, "src", "scRNAseq")
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "result", "scRNAseq", "figure")
infodir <- file.path(workdir, "result", "scRNAseq", "info")
df.infodir <- file.path(infodir, 'doubletfinder')
df.figdir <- file.path(figdir, 'doubletfinder')

# set a random seed
set.seed(98)

# load functions
source(file.path(srcdir, "my_functions.R"))

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# Load Seurat object
panc <- readRDS(file.path(infodir, 'panc.initial.rds'))

# doublet rate
# https://assets.ctfassets.net/an68im79xiti/1eX2FPdpeCgnCJtw4fj9Hx/7cb84edaa9eca04b607f9193162994de/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf
# according to 10X Genomics document, the doublet rate is correlated with the number of targeted/capture cells
# 1000 recovered cells ~ 0.8%
doublet.rate <- 0.008

# number of PCs
pcs <- 15

# resolution
res <- 0.8

# pN: number of generated artificial doublets, as a proportion of the merged real-artificial data
pN <- 0.25

# detect doublet per sample
for (sid in unique(panc$orig.ident)){
  print(paste("processing",sid))
  # subset Seurat object
  seu <- subset(panc, subset=orig.ident %in% c(sid))
  # pre-process Seurat object
  seu %<>% NormalizeData()
  seu %<>% FindVariableFeatures(selection.method="vst", nfeatures=2000)
  seu %<>% ScaleData()
  seu %<>% RunPCA()
  seu %<>% RunUMAP(dims=1:pcs)
  seu %<>% FindNeighbors(reduction="pca", dims=1:pcs)
  seu %<>% FindClusters(resolution=res)
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(seu, PCs=1:pcs, sct=FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  # factor --> numeric
  bcmvn$pK <- as.numeric(as.character(bcmvn$pK))
  write.table(bcmvn, file=file.path(df.infodir,paste("BCmvn","pK",sid,"txt",sep='.')), quote=F, sep='\t', row.names=F, col.names=T)
  # select optimal pK
  pK <- bcmvn$pK[which(bcmvn$BCmetric %in% max(bcmvn$BCmetric))]
  # Doublet Proportion Estimate
  nExp_poi <- round(doublet.rate*ncol(seu)/1000*ncol(seu))
  # run doublet detection
  seu %<>% doubletFinder_v3(PCs=1:pcs, pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=FALSE)
  # save pNN, singlet/doublet classification (original+adjusted), seurat_clusters (RNA_snn_res.0.8), UMAP coordinates(new)
  tdata <- FetchData(seu, vars=c(paste("pANN",pN,pK,nExp_poi,sep="_"),
                                 paste("DF.classifications",pN,pK,nExp_poi,sep="_"),
                                 'seurat_clusters','UMAP_1','UMAP_2'))
  write.table(tdata, file=file.path(df.infodir,paste("doublet","results",sid,"txt",sep='.')), quote=F, sep='\t', row.names=T, col.names=NA)
}

# collect doublet predictions
df.results <- data.frame(orig.ident=character(), Name=character(), Type=character(), Virus=character(), df.isdoublet=character())
for (sid in unique(panc.initial$orig.ident)){
  print(paste("processing",sid))
  # get cell info
  data <- FetchData(panc.initial, vars=c('orig.ident','Name','Type','Virus'), cells=grep(paste0('^',sid), rownames(panc.initial@meta.data), value=T))
  # read in DoubletFinder results
  tdf <- read.table(file.path(df.infodir, paste('doublet','results',sid,'txt',sep='.')), header=T, check.names=F, stringsAsFactors=F, sep='\t', row.names=1)
  colnames(tdf) <- c('df.score','df.isdoublet','seurat_clusters','UMAP_1','UMAP_2')
  # merge table
  data <- data %>% rownames_to_column('cellID') %>% 
    left_join(tdf %>% rownames_to_column('cellID') %>% select('cellID','df.isdoublet'), by='cellID') %>% 
    column_to_rownames('cellID')
  # add to results table
  df.results <- rbind(df.results, data)
}

# add DoubletFinder doublets labeling to seurat object
panc.initial %<>% AddMetaData(metadata=df.results[rownames(panc.initial@meta.data), c('df.isdoublet'),drop=F])

# save Seurat object (initial)
saveRDS(panc.initial, file=file.path(infodir, "panc.initial.rds"))
