# step_4.sample_integration.R
# step 4: integrate samples with MNN based correction.
# 

library(scran)
library(Seurat)
library(dplyr)
library(tibble)
library(magrittr)
library(future)
library(batchelor)
library(scater)

workdir <- "."
srcdir <- file.path(workdir, "src", "scRNAseq")
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "result", "scRNAseq", "figure")
infodir <- file.path(workdir, "result", "scRNAseq", "info")

# set a random seed
set.seed(98)

# load functions
source(file.path(srcdir, "my_functions.R"))

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# read in dissociation-related genes
dissofile <- paste(sourcedir, "dissociation", "dissociation_related_genes.human.txt", sep="/")
disso.genes <- as.vector(read.table(dissofile, header=F, check.names=F, sep="\t")$V1)

# Load Seurat object
panc.initial <- readRDS(file.path(infodir, 'panc.initial.rds'))

# get singlet cells (by DoubletFinder)
singlet.cells <- rownames(subset(FetchData(panc.initial, vars=c('df.isdoublet')), df.isdoublet == 'Singlet'))

# prepare raw UMI counts table from each donor
selected.donors <- c("XM1", "XM2", "XM3", "XM4")
sample.list <- list()
for (donor in selected.donors){
  sample.list[[donor]] <- panc.initial[["RNA"]]@counts[, intersect(rownames(subset(panc.initial@meta.data, orig.ident == donor)), singlet.cells)]
}

# create SingleCellExperiment object
sce.list <- list()
for (donor in names(sample.list)){
  sce.list[[donor]] <- SingleCellExperiment(list(counts=as.matrix(sample.list[[donor]])))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, min.size=200, assay.type="counts", method="hclust", min.mean=0.1))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) {computeSumFactors(x=x, min.mean=0.1, cluster=y)}, x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) {normalize(object=x)})

# rescale among donors
rescaled.sce.list <- do.call(multiBatchNorm, sce.list)

# create a seurat object with raw UMI counts
panc <- CreateSeuratObject(counts=as(do.call(cbind, lapply(rescaled.sce.list, function(x) counts(x))), "dgCMatrix"), 
                           project=project, assay="RNA", min.cells=0, min.features=0,
                           names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc[["percent.mito"]] <- PercentageFeatureSet(panc, pattern=mito.pattern)
panc[["percent.ribo"]] <- PercentageFeatureSet(panc, pattern=ribo.pattern)
panc[["percent.virus"]] <- PercentageFeatureSet(panc, pattern=virus.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc@meta.data[,"orig.ident"])])
}
panc %<>% AddMetaData(metadata=tmeta)

# replace normalized data with the scran normalized data
panc[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# Identification of highly variable features (feature selection)
panc %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc), c(disso.genes, 
                                                    grep(mito.pattern, rownames(panc), value=T), 
                                                    grep(ribo.pattern, rownames(panc), value=T), 
                                                    grep(virus.pattern, rownames(panc), value=T)))
variable.genes <- head(variable.genes, 3000)

# perform MMN-based correction
original <- lapply(rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc@meta.data),], key="MNN_", assay=DefaultAssay(panc))

# save Seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))

# free space
rm(preclust.list)
rm(panc.initial)
rm(sce.list)
rm(tmeta)
rm(tdic)
rm(mnn.out)
rm(original)
rm(rescaled.sce.list)
gc()
