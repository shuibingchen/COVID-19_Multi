# step_6.generate_plots.R
# step 6: generate plots in the manuscript
# 

library(Seurat)
library(dplyr)
library(tibble)
library(magrittr)
library(future)

workdir <- "."
srcdir <- file.path(workdir, "R", "scRNAseq")
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "result", "figure")
infodir <- file.path(workdir, "result", "info")

# set a random seed
set.seed(98)

# pattern for defining virus genes
virus.pattern <- "^CoV2"

# load functions
source(file.path(srcdir, "my_functions.R"))

# Load Seurat object
panc <- readRDS(file.path(infodir, 'panc.rds'))

# set color
my.cluster.color <- c('#fdbf6f','#1f78b4','#33a02c','#cab2d6','#e31a1c')
names(my.cluster.color) <- 0:4

myExpLowColor <- '#d9d9d9'
myExpHighColor <- '#b30000'

# sample info
sample.info <- data.frame(SeqName=c("CIART_Mock","CIART_SARS_CoV_2","WT_Mock","WT_SARS_CoV_2"), 
                          Name=c("CIART_Mock","CIART_SARS-CoV-2","WT_Mock","WT_SARS-CoV-2"), 
                          Type=c('CIART-/-','CIART-/-','WT','WT'), 
                          Virus=c('Mock','SARS-CoV-2','Mock','SARS-CoV-2'))
rownames(sample.info) <- c("XM1","XM2","XM3","XM4")

# Figure 2l: bar plot showing the percentage of virus UMI counts (relative to the total UMI counts) in each sample
stat.umi <- data.frame(virus.umi=colSums(panc[['RNA']]@counts[grep(virus.pattern, rownames(panc), value=T),]),
                       total.umi=colSums(panc[['RNA']]@counts))
stat.umi$sample <- sapply(rownames(stat.umi), function(x) { strsplit(x, '_')[[1]][1] })
stat.umi <- stat.umi %>% group_by(sample) %>% summarise_all(sum) %>% mutate(percent.virus.umi=round(virus.umi/total.umi*100,1)) %>%
  left_join(sample.info %>% rownames_to_column('sample'), by='sample')

stat.umi$Name <- factor(stat.umi$Name, levels=c("WT_Mock","CIART_Mock","WT_SARS-CoV-2","CIART_SARS-CoV-2"))
g <- ggplot(stat.umi, aes(x=Name, y=percent.virus.umi))
g <- g + geom_bar(stat="identity")
g <- g + ylab('Percentage of virus UMI counts')
g <- g + theme_bw()
g <- g + theme(axis.title.x = element_blank(), axis.title.y = element_text(size=15, face="bold", color='black'))
g <- g + theme(axis.text.x = element_text(size=14, face="bold", angle=60, hjust=1, color='black'), axis.text.y = element_text(size=14, face="bold", color='black'))
g <- g + theme(legend.position = "None")
ggsave(file.path(figdir, 'fig.2l.png'), plot=g, width=6, height=5.5, dpi=300)

# Figure 2m: UMAP illustrating the five cell clusters
g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", tcolor=my.cluster.color, tlabel=FALSE, tsplit=FALSE, tptsize=0.6) 
ggsave(file.path(figdir, "fig.2m.png"), plot=g, width=8.5, height=6, dpi=300)

# Load marker genes for each cluster
markers.wilcox.all <- read.table(file.path(infodir, "markers.pos.final_clusters.wilcox.minpct_0.25.txt"), header=T, sep='\t', row.names=1)

our.markers <- list()
for (k in 0:4){
  print(paste("load","our","cluster",k))
  # filter based on p_val_adj and avg_logFC
  tmarker <- subset(markers.wilcox.all, cluster==k & p_val_adj < 0.01 & avg_logFC > 0)
  our.markers[[paste0('Cluster_',k)]] <- rownames(tmarker)
}

# Load marker genes for the reference human adult lung dataset
ref.data.file <- file.path(sourcedir,'reference','ref.genes.txt.gz')
gene.list <- read.table(ref.data.file, header=T, check.names=F, stringsAsFactors=F, sep='\t')
ref.markers <- list()
for (c in unique(gene.list$Cluster)){
  ref.markers[[paste0('Ref_',c)]] <- as.vector(subset(gene.list, Cluster == c)[,'Gene'])
}

# Figure 2n: Correlation analysis of cell fates in our clusters and human adult lung cell clusters
plotFracOverlap(ref.markers, our.markers, infodir, figdir, "fig.2n", tsvg=T, twidth=7.5, theight=3.5)

# Figure 2o: percentage of virus UMIs per cell in WT and CIART-/- samples
data <- FetchData(panc, vars=c('ident','orig.ident','percent.virus')) %>% 
  filter(orig.ident %in% c('XM2','XM4')) %>% mutate(group=paste0('C',ident,'-',orig.ident))
# add a place holder for the missing state (no cells in cluster 2 for WT_SARS-CoV-2)
data <- rbind(data, data.frame(ident='2', orig.ident='XM4', percent.virus=NA, group='C2-XM4'))

# define orders in the plot
sample.order <- c()
for (k in 0:4){
  sample.order <- c(paste(paste0('C',k), c('XM2','XM4'), sep='-'), sample.order)
}
data$group <- factor(data$group, levels=sample.order)

g <- ggplot(data, aes(x=group, y=percent.virus))
g <- g + geom_violin(aes(fill=orig.ident), scale="width", alpha=0.7)
g <- g + geom_jitter(aes(color=orig.ident), shape=16, position=position_jitter(width=0.2, height=0), alpha=0.2)
g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray20")
g <- g + scale_color_manual(values=c('XM2'='#a6cee3','XM4'="#fb9a99"))
g <- g + scale_fill_manual(values=c('XM2'='#1f78b4','XM4'='#e31a1c'))
g <- g + coord_flip(ylim=c(0,100))
g <- g + ylab('Percentage of virus UMIs')
g <- g + theme_bw()
g <- g + theme(axis.text.x=element_text(color='black'), axis.text.y=element_text(color='black'))
g <- g + theme(axis.title.x=element_text(color='black'), axis.title.y=element_blank())
g <- g + theme(legend.position="none")
ggsave(file.path(figdir, 'fig.2p.png'), plot=g, width=3, height=6.5, dpi=300)

# Figure 2p: dot plot showing the expression of all SARS-CoV-2 viral genes in WT and CIART-/- infected samples in the "Ciliated-like cells 1" cluster
g <- DotPlot(subset(panc, subset=final.clust %in% c(1) & orig.ident %in% c('XM2','XM4')), assay='RNA',
             features=grep(virus.pattern, rownames(panc), value=T), group.by='Name')
g <- g + scale_color_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B")
g <- g + coord_flip()
g <- g + theme_bw()
g <- g + theme(axis.text.x=element_text(angle=60, hjust=1, color='black'), axis.text.y=element_text(color='black'))
g <- g + theme(axis.title=element_blank())
ggsave(file.path(figdir, "fig.2p.png"), plot=g, width=4, height=6, dpi=300)

# Extended Data Figure 6a: UMAP plot showing the expression of selected marker genes
markers <- c('FOXJ1','KRT5','TP63','MKI67','TOP2A','CHGA')
for (gene in markers){
  print(gene)
  plot <- MyFeaturePlot(tobj=panc, tgenes=gene, tcells=NULL, tassay="RNA", treduction.name="umap", tsort=TRUE, tncol=1, 
                        tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tlegend=NULL)
  ggsave(file.path(figdir, paste('extended.fig.6a',gene,'png',sep='.')), plot=plot, width=8, height=6, dpi=300)
}

# Extended Data Figure 6b: stacked bar plot showing cell compositions in each sample
g <- myStackBarCellComposition(panc, tcluster_by='ident', tcluster_order=0:4, tgroup_by='Name',
                               tgroup_order=c('WT_Mock','WT_SARS-CoV-2','CIART_Mock','CIART_SARS-CoV-2'),
                               tcells=NULL, tanncolor=my.cluster.color, ttlsize=20, ttxsize=18, tltsize=16)
ggsave(file.path(figdir, 'extended.fig.6b.png'), plot=g, width=6, height=6.5)
