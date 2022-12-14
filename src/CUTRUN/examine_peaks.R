# examine_peaks.R
# review the peaks in CUT&RUN samples
# 

library(DiffBind, quietly=T, warn.conflicts=FALSE)
library(ChIPseeker, quietly=T, warn.conflicts=FALSE)
library(ggplot2, quietly=T, warn.conflicts=FALSE)
library(dplyr, quietly=T, warn.conflicts=FALSE)
library(tidyr, quietly=T, warn.conflicts=FALSE)
library(tibble, quietly=T, warn.conflicts=FALSE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly=T, warn.conflicts=FALSE)
library(org.Hs.eg.db, quietly=T, warn.conflicts=FALSE)

workdir <- '.'
sourcedir <- file.path(workdir, "data")

# narrowPeak files can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202966
peakdir <- file.path(sourcedir, "cutrun", "peaks")

figdir <- file.path(workdir, "results", "cutrun", "figure")
infodir <- file.path(workdir, "results", "cutrun", "info")

# create a sample sheet
sample.ids <- c("WT_1", "WT_2", "IgG_Control")
samples <- data.frame(SampleID=sample.ids,
    Condition=c(rep('WT',2), 'IgG'),
    Replicate=c(1,2,1),
    Peaks=file.path(peakdir, paste(sample.ids, 'narrowPeak', sep='.')),
    PeakCaller=rep('narrow',3))

# create a dba object
panc <- dba(minOverlap=1, sampleSheet=samples)

# filter peaks:
# - For WT condition, take the intersecting peaks from the two replicates
# - Mask any peaks observed in the IgG control sample
panc_consensus <- dba.peakset(panc, consensus=DBA_CONDITION, minOverlap=2)
panc_OL <- dba.overlap(panc_consensus, panc_consensus$masks$Consensus | panc_consensus$masks$IgG)

# annotate peaks
# load the knownGene track
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# prepare the TSS regions, which are defined as the flanking sequence of the TSS sites.
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# prepare peak sets
peak.list <- list(candidate=panc_OL$onlyB)

# collect peak data
peakAnnoList <- lapply(peak.list, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), annoDb="org.Hs.eg.db", verbose=FALSE)

# align the peak that are mapping to the TSS regions, and generate the tagMatrix
tagMatrixList <- lapply(peak.list, getTagMatrix, windows=promoter)

# Pie chart to visualize the genomic annotation
for (k in 1:length(peakAnnoList)){
  myset <- names(peakAnnoList)[k]
  print(myset)
  peakAnno = peakAnnoList[[k]]
  png(file.path(figdir, 'piechart.peaks.png'), width = 6400, height = 4000, units = "px", pointsize = 6, res=600)
  print(plotAnnoPie(peakAnno))
  dev.off()

  # convert results to data frame
  myresult <- as.data.frame(peakAnno)
  # add a column of simplified annotation
  myresult$simple.annotation <- sapply(myresult$annotation, 
                                       function(x) { strsplit(x, '[ (]')[[1]][1] })
  # write to file
  write.table(myresult, 
              file.path(infodir, paste('annotated','peaks',myset,'txt',sep='.')), 
              quote=F, sep='\t', row.names=F)
}

# Average Profile of ChIP peaks binding to TSS region
png(file.path(figdir, 'average.profile.peaks.png'), width = 6400, height = 4000, units = "px", pointsize = 6, res=600)
plotAvgProf(tagMatrixList[[1]], xlim=c(-3000, 3000), color='red') + theme(legend.position="none")
dev.off()

# Peak heatmap of ChIp peaks binding to TSS region
png(file.path(figdir, 'heatmap.peaks.png'), width = 6400, height = 4000, units = "px", pointsize = 6, res=600)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), title=NULL, color='red')
dev.off()

sessionInfo()

