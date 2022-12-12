# step_5.cell_clustering.R
# step 5: cluster cells and infer cell types
# 

library(Seurat)
library(dplyr)
library(tibble)
library(magrittr)
library(future)

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

# Load Seurat object
panc <- readRDS(file.path(infodir, 'panc.rds'))

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc %<>% FindNeighbors(reduction="mnn", dims=1:50)
# FindClusters
panc %<>% FindClusters(resolution=seq(0.05,1,by=0.05), verbose=T)

# Run non-linear dimensional reduction (UMAP)
panc %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=35, n.epochs=2500)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.2")

# reorder clusters
Idents(panc) <- factor(Idents(panc), levels=0:(length(unique(Idents(panc)))-1))

# Find markers for each cluster
markers.wilcox.all <- FindAllMarkers(panc, only.pos=TRUE, logfc.threshold=0.25, test.use="wilcox", min.pct=0.25, assay="RNA")

# save current cluster
panc[["base.clust"]] <- Idents(panc)

# merge clusters after manual review
# C1 + C3 + C5 + C6 + C8 + C11 + C12 + C14 + C16            ==>  C0 (Basal cells)
# C2 + C4 + C10 + C13                                       ==>  C1 (Ciliated-like cells - 1)
# C0                                                        ==>  C2 (Ciliated-like cells - 2)
# C7 + C9                                                   ==>  C3 (Proliferating basal cells)
# C15                                                       ==>  C4 (Neuronendocrine)
merge.clust <- c(2,0,1,0,1,0,0,3,0,3,1,0,0,1,0,4,0)
names(merge.clust) <- 0:16

# perform hierarchical clustering on the initial clusters based on marker gene similarity to the reference human adult lung data
# collect marker genes of our clusters
our.markers <- list()
for (k in 0:16){
  print(paste("load","our","cluster",k))
  # filter based on p_val_adj and avg_logFC
  tmarker <- subset(markers.wilcox.all, cluster==k & p_val_adj < 0.01 & avg_logFC > 0)
  our.markers[[paste0('Our_',k)]] <- rownames(tmarker)
}

# load marker genes of the reference set
ref.data.file <- file.path(sourcedir,'reference','ref.genes.txt.gz')
gene.list <- read.table(ref.data.file, header=T, check.names=F, stringsAsFactors=F, sep='\t')
ref.markers <- list()
for (c in unique(gene.list$Cluster)){
  ref.markers[[paste0('Ref_',c)]] <- as.vector(subset(gene.list, Cluster == c)[,'Gene'])
}

# prepare annotations in the heatmap
annrows <- data.frame(merged.clust=merge.clust)
rownames(annrows) <- paste('Cluster',0:16,sep='_')

# set color
my.cluster.color <- c('#fdbf6f','#1f78b4','#33a02c','#cab2d6','#e31a1c')
names(my.cluster.color) <- 0:4
anncolors <- list('merged.clust'=my.cluster.color)

# generate heatmap based on marker genes similarity and perform hierarchical clustering on the initial clusters
plotFracOverlap(ref.markers, our.markers, infodir, figdir, "heatmap.marker.similarity", 
                tsvg=F, cluster_rows=T, clustering_method="complete", annrows=annrows, anncolors=anncolors, twidth=5.5, theight=4)

# merge clusters based on the hierarchical clustering result
final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:4)

# add final clusters to meta data
panc[["final.clust"]] <- final.clust

# set final clusters
Idents(panc) <- final.clust

# Find markers for each cluster
markers.wilcox.all <- FindAllMarkers(panc, only.pos=TRUE, logfc.threshold=0.25, test.use="wilcox", min.pct=0.25, assay="RNA")
write.table(as.data.frame(markers.wilcox.all), file=file.path(infodir, "markers.pos.final_clusters.wilcox.minpct_0.25.txt"), quote=FALSE, na="", sep="\t", col.names=NA)

# save seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))
