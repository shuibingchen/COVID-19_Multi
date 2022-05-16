# run_cluster.R
# perform clustering analysis on scRNAseq data
# Author: Tuo Zhang
# Date: 03/10/2021
# Version: 1.0
# 

library(scran)
library(Seurat)
library(dplyr)
library(tibble)
library(magrittr)
library(MAST)
library(future)
library(batchelor)
library(scater)
library(ggsci)
library(RColorBrewer)

workdir <- "."
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

# project name
project <- "shuibing"

# pattern for defining mitochondrial/ribosomal/virus genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"
virus.pattern <- "^CoV2"

# load functions
setwd(workdir)
source("my_functions.R")

# set a random seed
set.seed(98)

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# read in dissociation-related genes
dissofile <- paste(sourcedir, "dissociation", "dissociation_related_genes.human.txt", sep="/")
disso.genes <- as.vector(read.table(dissofile, header=F, check.names=F, sep="\t")$V1)

# ----------------------------------------------------- filter cells ---------------------------------------------------- #
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
  raw.counts.list[[k]] <- my.Read10X(file.path(sourcedir, sid, "filtered_feature_bc_matrix"), pid)
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

# perform cell filtering
# nGene > 500, nGene <= 8000, nUMI > 1000, nUMI <= 70000, percent.mito < 15%
panc.initial %<>% subset(subset=nFeature_RNA_nonVirus > 500 & nFeature_RNA_nonVirus <= 8000 & nCount_RNA_nonVirus > 1000 & nCount_RNA_nonVirus <= 70000 & percent.mito < 15)

# free space
rm(raw.counts.list)
rm(raw.counts.all)
rm(tmeta)
rm(tdic)
# ----------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------- run MNN-based correction ----------------------------------------------- #
# prepare raw UMI counts table from each donor
selected.donors <- c("XM1", "XM2", "XM3", "XM4")
sample.list <- list()
for (donor in selected.donors){
  sample.list[[donor]] <- panc.initial[["RNA"]]@counts[, rownames(subset(panc.initial@meta.data, orig.ident == donor))]
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

# remove dissociation-related genes, mitochondrial/ribosomal/virus genes from variable gene list
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

# free space
rm(panc.initial)
rm(preclust.list)
rm(sce.list)
rm(rescaled.sce.list)
rm(original)
rm(mnn.out)
rm(tmeta)
rm(tdic)
# ----------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------- run UMAP and clustering ----------------------------------------------- #
# Run non-linear dimensional reduction (UMAP)
panc %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=35, n.epochs=2500)

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc %<>% FindNeighbors(reduction="mnn", dims=1:50)
# FindClusters
panc %<>% FindClusters(resolution=0.45, verbose=T)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.45")

# reorder clusters
Idents(panc) <- factor(Idents(panc), levels=0:(length(unique(Idents(panc)))-1))

# Finding differentially expressed features (cluster biomarkers)
markers.wilcox.all <- FindAllMarkers(panc, only.pos=TRUE, logfc.threshold=0.25, test.use="wilcox", min.pct=0.25, assay="RNA")
# ----------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------------- merge clusters ---------------------------------------------------- #
# save current cluster
panc[["base.clust"]] <- Idents(panc)

# merge clusters after manual review
# C0 + C1 + C2 + C5 + C6 + C8 + C10 + C11 + C12 + C13 + C14 + C17 + C18 + C20   ==>  C0 (Ciliated like cells)
# C3 + C4 + C16                                                                 ==>  C1 (Basal cells)
# C7 + C9                                                                       ==>  C2 (Proliferating cells)
# C15                                                                           ==>  C3 (Unknown cells)
# C19                                                                           ==>  C4 (Neuroendocrine cells)
merge.clust <- c(0,0,0,1,1,0,0,2,0,2,0,0,0,0,0,3,1,0,0,4,0)
names(merge.clust) <- 0:20

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:4)

# add final clusters to meta data
panc[["final.clust"]] <- final.clust

# set final clusters
Idents(panc) <- final.clust

# Finding differentially expressed features (cluster biomarkers)
for (k in c(0:4)){
  print(paste("Cluster",k))
  # wilcox test
  myDETest(panc, k,"wilcox",TRUE,infodir,figdir,tassay="RNA",tsuf="merged_clusters")
}

# save seurat object
saveRDS(panc, file=paste(infodir, "panc.rds", sep="/"))

# free memory
rm(merge.clust)
rm(final.clust)
# ----------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------- generate figures ------------------------------------------------ #
# Figure 2l: bar plot showing the percentage of virus UMI counts (relative to the total UMI counts) in each sample
stat.umi <- data.frame(virus.umi=colSums(panc[['RNA']]@counts[grep(virus.pattern, rownames(panc), value=T),]),
                       total.umi=colSums(panc[['RNA']]@counts))
stat.umi$sample <- sapply(rownames(umis), function(x) { strsplit(x, '_')[[1]][1] })
stat.umi <- stat.umi %>% group_by(sample) %>% summarise_all(sum) %>% mutate(percent.virus.umi=round(virus.umi/total.umi*100,1)) %>%
  left_join(sample.info %>% rownames_to_column('sample'), by='sample')

stat.umi$Name <- factor(stat.umi$Name, levels=c("WT_Mock","WT_SARS-CoV-2","CIART_Mock","CIART_SARS-CoV-2"))
g <- ggplot(stat.umi, aes(x=Name, y=percent.virus.umi))
g <- g + geom_bar(stat="identity")
g <- g + ylab('Percentage of virus UMI counts')
g <- g + theme_bw()
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_text(size=15, face="bold", color='black'))
g <- g + theme(axis.text.x = element_text(size=14, face="bold", angle=60, hjust=1, color='black'), axis.text.y = element_text(size=14, face="bold", color='black'))
g <- g + theme(legend.position = "None")
ggsave(file.path(figdir,'percent.virus.umi.per.sample.png'), plot=g, width=6, height=5.5, dpi=300)

# Figure 2m: UMAP illustrating the five cell clusters
my.cluster.color <- c('#fdbf6f','#1f78b4','#cab2d6','#33a02c','#e31a1c')
names(my.cluster.color) <- 0:4
g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", tcolor=my.cluster.color, tlabel=FALSE, tsplit=FALSE, tptsize=0.6) 
ggsave(file.path(figdir, "UMAPPlot.by_Cluster.png"), plot=g, width=8.5, height=6, dpi=300)

#######################################
# compare to known reference data set #
#######################################
# load the reference marker genes
ref.data.file <- file.path(sourcedir,'reference','ref.genes.txt.gz')
gene.list <- read.table(ref.data.file, header=T, check.names=F, stringsAsFactors=F, sep='\t')
ref.markers <- list()
for (c in unique(gene.list$Cluster)){
  ref.markers[[paste0('Ref_',c)]] <- as.vector(subset(gene.list, Cluster == c)[,'Gene'])
}

# load the marker genes for our clusters
our.markers <- list()
# cutoff for p_val_adj
padj.cut <- 0.01
for (k in 0:4){
  print(paste("load","our","cluster",k))
  # load marker
  tmarker <- read.table(file.path(infodir, paste0("C",k,".wilcox.merged_clusters.origin.markers.txt")), header=T, check.names=F, stringsAsFactors=F, row.names=1, sep='\t')
  # filter based on p_val_adj and avg_logFC
  tmarker <- subset(tmarker, p_val_adj < padj.cut & avg_logFC > 0)
  our.markers[[paste0('Our_',k)]] <- rownames(tmarker)
}

# Figure 2n: Enrichement analysis of our dataset using genes highly expressed in adult human ciliates or proximal ciliated cells
scoreByRefMarkers(panc, ref.markers, trefPattern='^Ref_', infodir, figdir, 'enrichment.ref_cluster', tsvg=T)

# Figure 2o: Correlation analysis of genes with cell fates in our dataset and adult human lung cells
plotFracOverlap(ref.markers, our.markers, infodir, figdir, "correlation.ref_cluster", tsvg=T)

# Figure 2p: percentage of virus UMIs in WT and CIART-/- samples
data <- FetchData(panc, vars=c('ident','orig.ident','percent.virus')) %>% 
  filter(orig.ident %in% c('XM2','XM4'))
# reorder
sample.order <- c()
for (k in 0:4){
  sample.order <- c(paste(paste0('C',k), c('XM2','XM4'), sep='-'), sample.order)
}
data$ident <- factor(data$ident, levels=sample.order)

g <- ggplot(data, aes(x=ident, y=percent.virus))
g <- g + geom_violin(aes(fill=orig.ident), scale="width", alpha=0.7)
g <- g + geom_jitter(aes(color=orig.ident), shape=16, position=position_jitter(width=0.2, height=0), alpha=0.2)
g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray20")
g <- g + scale_color_manual(values=c('XM2'='#a6cee3','XM4'="#fb9a99"))
g <- g + scale_fill_manual(values=c('XM2'='#1f78b4','XM4'='#e31a1c'))
g <- g + coord_flip()
g <- g + ylab('Percentage of virus UMIs')
g <- g + theme_bw()
g <- g + theme(axis.text.x=element_text(color='black'), axis.text.y=element_text(color='black'))
g <- g + theme(axis.title.x=element_text(color='black'), axis.title.y=element_blank())
g <- g + theme(legend.position="none")
ggsave(file.path(figdir, 'virus.content.all_clusters.KO_SARS-CoV-2.vs.WT_SARS-CoV-2.png'), plot=g, width=3, height=5.5, dpi=300)

# Figure 2q: percent of infected cells in WT and CIART-/- samples
data <- FetchData(panc, vars=c('ident','orig.ident','percent.virus')) %>% 
  filter(orig.ident %in% c('XM2','XM4')) %>% 
  group_by(ident) %>% summarise_at('percent.virus', function(x) { round(sum(x>20)/length(x)*100,2) }) %>% 
  left_join(FetchData(panc, vars=c('ident','orig.ident')) %>% distinct(), by='ident') %>% 
  rename('percent.virus'='percent.infect')
data$ident <- factor(data$ident, levels=sample.order)

g <- ggplot(data, aes(x=ident, y=percent.infect))
g <- g + geom_bar(aes(fill=orig.ident), stat='identity')
g <- g + scale_fill_manual(values=c('XM2'='#1f78b4','XM4'='#e31a1c'))
g <- g + coord_flip()
g <- g + ylab('Percentage of infected cells')
g <- g + theme_bw()
g <- g + theme(axis.text.x=element_text(color='black'), axis.text.y=element_text(color='black'))
g <- g + theme(axis.title.x=element_text(color='black'), axis.title.y=element_blank())
g <- g + theme(legend.position="none")
ggsave(file.path(figdir, 'percent.infected.cells.cut_20.all_clusters.KO_SARS-CoV-2.vs.WT_SARS-CoV-2.png'), plot=g, width=3, height=5.5, dpi=300)

# Figure 2r: Dot plot showing the expression of SARS-CoV-2 viral genes in the Ciliated-like cell cluster (0) in the WT and CIART-/- infected samples 
g <- DotPlot(subset(panc, subset=final.clust %in% c(0) & orig.ident %in% c('XM2','XM4')), assay='RNA', 
             features=grep(virus.pattern, rownames(panc), value=T), group.by='Name')
g <- g + scale_color_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B")
g <- g + coord_flip()
g <- g + theme_bw()
g <- g + theme(axis.text.x=element_text(angle=60, hjust=1, color='black'), axis.text.y=element_text(color='black'))
g <- g + theme(axis.title=element_blank())
ggsave(file.path(figdir, 'exp.Dot.covid.genes.C0.KO_SARS-CoV-2.vs.WT_SARS-CoV-2.png'), plot=g, width=4, height=6, dpi=300)

# Extended Data Figure 6a: UMAP plots showing the expression in several marker genes
markers <- c('FOXJ1','TP63','KRT5','MKI67','TOP2A','CHGA')

myExpLowColor <- '#d9d9d9'
myExpHighColor <- '#b30000'

for (gene in markers){
  print(gene)
  plot <- MyFeaturePlot(tobj=panc, tgenes=gene, tcells=NULL, tassay="RNA", treduction.name="umap", tsort=TRUE, tncol=1, tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tlegend=NULL, twidth=15, theight=12.5, tunits="in", tres=300)
  ggsave(file.path(figdir, paste('exp','UMAP',gene,'png',sep='.')), plot=plot, width=8, height=6, dpi=300)
}

# Extended Data Figure 6b: Bar plot showing the composition of the five cell clusters in each sample
g <- myStackBarCellComposition(panc, tcluster_by='ident', tcluster_order=0:4, tgroup_by='Name', 
                               tgroup_order=c("WT_Mock","WT_SARS-CoV-2","CIART_Mock","CIART_SARS-CoV-2"), 
                               tcells=NULL, tanncolor=my.cluster.color, ttlsize=20, ttxsize=18, tltsize=16)
ggsave(file.path(refined.figdir.2, 'percentage.cells.per_cluster.by_sample.png'), plot=g, width=6, height=6.5, dpi=300)
# ----------------------------------------------------------------------------------------------------------------------- #
