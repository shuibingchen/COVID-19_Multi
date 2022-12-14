# 1.identify_DE_genes_in_SARS-CoV-2.R
# perform differential expression (DE) analysis to identify genes involved in SARS-CoV-2 infection in multiple organs
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
perform_de_analysis <- function(countsfile, cdt.order, alpha, infodir, figdir, suffix, ann.file, num.cores=4){
    # read in raw counts file
    countData <- read.table(file=countsfile, header=TRUE, check.names=FALSE, stringsAsFactors=F, sep='\t', row.names=1)
    # use human genes only
    countData <- countData[grep("^GU280_",rownames(countData),value=T,invert=T),]
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
    for (i in c(2:length(cdt.order))){
        print(paste(cdt.order[i], cdt.order[1], sep=" v.s. "))
        res <- results(dds, contrast=c("condition", cdt.order[i], cdt.order[1]), alpha=alpha, parallel=TRUE, BPPARAM=MulticoreParam(num.cores))
        # moderated log2 fold changes
        resLFC <- lfcShrink(dds, res=res, type="ashr", parallel=TRUE, BPPARAM=MulticoreParam(num.cores))
        # reorder results by adjusted-pvalue
        resLFCOrdered <- resLFC[order(resLFC$padj),]
        # add annotations
        resLFCOrdered <- merge(resLFCOrdered, ann, by='row.names', all.x=T)
        colnames(resLFCOrdered) <- c('gene', colnames(resLFCOrdered)[-1])
        # exporting results to file
        write.table(as.data.frame(resLFCOrdered), file=file.path(infodir, paste("DE", cdt.order[i], "vs", cdt.order[1], "LFC", "tsv", sep=".")), quote=FALSE, sep='\t', row.names=F)
    }
    # regularized log transform counts for visualization
    rld <- rlog(dds, blind=F)
    # add annotations
    rld.data <- as.data.frame(assay(rld))
    rld.data <- merge(rld.data, ann, by='row.names', all.x=T)
    colnames(rld.data) <- c('gene', colnames(rld.data)[-1])
    # save to file
    write.table(rld.data, file=file.path(infodir, paste("rld",suffix,"tsv", sep=".")), quote=FALSE, sep='\t', col.names=NA)
}

# false discovery rate cutoff
alpha <- 0.1

# use parallelization
register(MulticoreParam(4))

# counts files can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202963
# Airway 48h Infected vs. Mock
perform_de_analysis(countsfile=file.path(countsdir, 'counts.AWO.48h.txt'), 
    cdt.order=c('Airway.48h.Mock','Airway.48h.MOI0.01','Airway.48h.MOI0.1','Airway.48h.MOI1'), 
    alpha=alpha, infodir=infodir, figdir=figdir, suffix="Airway.48h", 
    ann.file=file.path(sourcedir, 'annotation', 'grch37.ensembl.genes.txt.gz'), num.cores=4)

# Alveolar 48h Infected vs. Mock
perform_de_analysis(countsfile=file.path(countsdir, 'counts.ALO.48h.txt'), 
    cdt.order=c('Alveolar.48h.Mock','Alveolar.48h.MOI0.01','Alveolar.48h.MOI0.1','Alveolar.48h.MOI1'), 
    alpha=alpha, infodir=infodir, figdir=figdir, suffix="Alveolar.48h", 
    ann.file=file.path(sourcedir, 'annotation', 'grch37.ensembl.genes.txt.gz'), num.cores=4)

# Cardiomyocyte 48h Infected vs. Mock
perform_de_analysis(countsfile=file.path(countsdir, 'counts.CM.48h.txt'), 
    cdt.order=c('Cardiomyocyte.48h.Mock','Cardiomyocyte.48h.MOI0.01','Cardiomyocyte.48h.MOI0.1','Cardiomyocyte.48h.MOI1'), 
    alpha=alpha, infodir=infodir, figdir=figdir, suffix="Cardiomyocyte.48h", 
    ann.file=file.path(sourcedir, 'annotation', 'grch37.ensembl.genes.txt.gz'), num.cores=4)

sessionInfo()

