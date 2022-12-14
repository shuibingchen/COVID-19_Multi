# 4.identify_DE_genes_in_CIART_KO.R
# perform differential expression (DE) analysis to identify genes involved in CIART-/- in multiple organs
# 

library(DESeq2)
library(BiocParallel)
library(vsn)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

workdir <- '.'
sourcedir <- file.path(workdir, "data")
countsdir <- file.path(sourcedir, "bulk_counts")
figdir <- file.path(workdir, "results", "bulkRNAseq", "figure")
infodir <- file.path(workdir, "results", "bulkRNAseq", "info")

# function to process data from one organ
perform_de_analysis <- function(countsfile, cdt.order, cmp.list, alpha, infodir, figdir, suffix, ann.file, num.cores=4){
    # read in raw counts file
    countData <- read.table(file=countsfile, header=TRUE, check.names=FALSE, stringsAsFactors=F, sep='\t', row.names=1)
    # create experimental label
    colData <- data.frame(condition=sapply(strsplit(colnames(countData), '_'), FUN=function(x) { paste(x[-length(x)], collapse='.') }))
    rownames(colData) <- colnames(countData)
    colData$condition <- factor(colData$condition, levels=cdt.order)
    # construct a DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~condition)
    # update factor levels
    dds$condition <- factor(dds$condition, levels=cdt.order)
    # run DESeq2
    dds <- DESeq(dds, parallel=TRUE, BPPARAM=MulticoreParam(num.cores))
    # load gene annotations
    ann <- read.table(ann.file, header=T, check.names=F, stringsAsFactors=F, sep='\t', row.names=1)
    # collect DE results
    for (cmp.pair in cmp.list){
        print(paste(cmp.pair[2], cmp.pair[1], sep=" v.s. "))
        res <- results(dds, contrast=c("condition", cmp.pair[2], cmp.pair[1]), alpha=alpha, parallel=TRUE, BPPARAM=MulticoreParam(num.cores))
        # moderated log2 fold changes
        resLFC <- lfcShrink(dds, res=res, type="ashr", parallel=TRUE, BPPARAM=MulticoreParam(num.cores))
        # reorder results by adjusted-pvalue
        resLFCOrdered <- resLFC[order(resLFC$padj),]
        # add annotations
        resLFCOrdered <- merge(resLFCOrdered, ann, by='row.names', all.x=T)
        colnames(resLFCOrdered) <- c('gene', colnames(resLFCOrdered)[-1])
        # exporting results to file
        write.table(as.data.frame(resLFCOrdered), file=file.path(infodir, paste("DE", cmp.pair[2], "vs", cmp.pair[1], "LFC", "tsv", sep=".")), quote=FALSE, sep='\t', row.names=F)
    }
    # rlog transform counts for visualization
    rld <- rlog(dds, blind=F)
    # add annotations
    rld.data <- as.data.frame(assay(rld))
    rld.data <- merge(rld.data, ann, by='row.names', all.x=T)
    colnames(rld.data) <- c('gene', colnames(rld.data)[-1])
    # save to file
    write.table(rld.data, file=file.path(infodir, paste("rld",suffix,"tsv", sep=".")), quote=FALSE, sep='\t', col.names=NA)
    # select human genes
    rld.noCov <- rld[grep("^GU280_",rownames(rld),value=T,invert=T),]
    for (cmp.pair in cmp.list){
        # use samples in the comparison pair
        tdata <- rld.noCov[,grep(paste(cmp.pair, collapse='|'), colnames(rld.noCov), value=T)]
        # perform PCA
        mypcadata <- plotPCA(tdata, intgroup=c('condition'), returnData=TRUE)
        percentVar <- round(100*attr(mypcadata, "percentVar"))
        g <- ggplot(mypcadata, aes(PC1, PC2, color=condition, label=name))
        g <- g + geom_point(size=3, shape=19) + xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance"))
        g <- g + theme_classic()
        ggsave(file=file.path(figdir, paste("PCA",suffix,"png",sep='.')), plot=g, width = 8, height = 6, dpi = 600)
        # heatmap of the sample-to-sample distances
        sampleDists <- dist(t(assay(tdata)))
        sampleDistMatrix <- as.matrix(sampleDists)
        hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(255)
        png(file=file.path(figdir, paste("heatmap.sample-to-sample",suffix,"png",sep='.')), width = 3600, height = 3200, units = "px", pointsize = 6, res=600)
        pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=hmcol, fontsize_row=6, fontsize_col=6)
        dev.off()
    }
}

# false discovery rate cutoff
alpha <- 0.1

# use parallelization
register(MulticoreParam(4))

# counts files can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202963
# Airway: CIART-/- vs. WT
perform_de_analysis(countsfile=file.path(countsdir, 'counts.AWO.downsample.txt'), 
    cdt.order=c('Airway.WT.Mock','Airway.WT.MOI0.1','Airway.CIART.Mock','Airway.CIART.MOI0.1'), 
    cmp.list=list(c('Airway.WT.Mock','Airway.CIART.Mock'), c('Airway.WT.MOI0.1','Airway.CIART.MOI0.1')), 
    alpha=alpha, infodir=infodir, figdir=figdir, suffix="Airway.CIART", 
    ann.file=file.path(sourcedir, 'annotation', 'grch37.ensembl.genes.txt.gz'), num.cores=4)

# Alveolar: CIART-/- vs. WT
perform_de_analysis(countsfile=file.path(countsdir, 'counts.ALO.downsample.txt'), 
    cdt.order=c('Alveolar.WT.Mock','Alveolar.WT.MOI0.1','Alveolar.CIART.Mock','Alveolar.CIART.MOI0.1'), 
    cmp.list=list(c('Alveolar.WT.Mock','Alveolar.CIART.Mock'), c('Alveolar.WT.MOI0.1','Alveolar.CIART.MOI0.1')), 
    alpha=alpha, infodir=infodir, figdir=figdir, suffix="Alveolar.CIART", 
    ann.file=file.path(sourcedir, 'annotation', 'grch37.ensembl.genes.txt.gz'), num.cores=4)

# Cardiomyocyte: CIART-/- vs. WT
perform_de_analysis(countsfile=file.path(countsdir, 'counts.CM.downsample.txt'), 
    cdt.order=c('Cardiomyocyte.WT.Mock','Cardiomyocyte.WT.MOI0.1','Cardiomyocyte.CIART.Mock','Cardiomyocyte.CIART.MOI0.1'), 
    cmp.list=list(c('Cardiomyocyte.WT.Mock','Cardiomyocyte.CIART.Mock'), c('Cardiomyocyte.WT.MOI0.1','Cardiomyocyte.CIART.MOI0.1')), 
    alpha=alpha, infodir=infodir, figdir=figdir, suffix="Cardiomyocyte.CIART", 
    ann.file=file.path(sourcedir, 'annotation', 'grch37.ensembl.genes.txt.gz'), num.cores=4)

sessionInfo()

