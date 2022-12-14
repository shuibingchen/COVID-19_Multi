# identify_differential_binding_sites.R
# perform differential binding analysis to identify differential binding sites in CIART-/- compared to WT
# 

library(DiffBind)
library(ggplot2, quietly=T, warn.conflicts=FALSE)
library(dplyr, quietly=T, warn.conflicts=FALSE)
library(tidyr, quietly=T, warn.conflicts=FALSE)
library(tibble, quietly=T, warn.conflicts=FALSE)
library(profileplyr, quietly=T, warn.conflicts=FALSE)
library(csaw, quietly=T, warn.conflicts=FALSE)
library(ChIPseeker, quietly=T, warn.conflicts=FALSE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, quietly=T, warn.conflicts=FALSE)
library(org.Hs.eg.db, quietly=T, warn.conflicts=FALSE)

workdir <- '.'
sourcedir <- file.path(workdir, "data")

# narrowPeak files can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202965
peakdir <- file.path(sourcedir, "atac", "peaks")
# bam files can be generated using our ATACseq pipeline: https://github.com/freshtuo/GRCF_ATACseq
bamdir <- file.path(sourcedir, "atac", "bam")

figdir <- file.path(workdir, "results", "ATACseq", "figure")
infodir <- file.path(workdir, "results", "ATACseq", "info")

# create a sample sheet
sample.ids <- c(paste('CM_WT',1:3,sep='_'), paste('CM_CIART',1:3,sep='_'))
samples <- data.frame(SampleID=sample.ids, 
    Condition=c(rep('WT',3), rep('CIART',3)),
    Replicate=rep(1:3,2), 
    bamReads=file.path(bamdir, paste(sample.ids, 'bam', sep='.')), 
    Peaks=file.path(peakdir, paste(sample.ids, 'narrowPeak', sep='.')),
    PeakCaller=rep('narrow',6))

# Perform a quick occupancy analysis on the peaksets
# Read in peaksets, generate a consensus peakset by: 
# - calculate consensus peaks overlap in at least two replicate samples for each condition 
# - calculate the final consensus peaks set by taking a union of peaks from all conditions 
# - Discard a pre-defined list of regions, as part of the ENCODE project, specific to the reference genome that are known to be problematic.

# create a dba object
panc <- dba(minOverlap=1, sampleSheet=samples)

# calculate consensus peaks overlap in at least two replicate samples in each condition
# this will add condition-level masks
panc_consensus <- dba.peakset(panc, consensus=DBA_CONDITION, minOverlap=2)

# calculate the final consensus peaks set
panc_consensus <- dba(panc_consensus, mask=panc_consensus$masks$Consensus, minOverlap=1)

# apply ENCODE blacklist
panc_consensus <- dba.blacklist(panc_consensus, blacklist=DBA_BLACKLIST_GRCH37, greylist=FALSE)

# extract consensus peaks for reads counting
consensus.peaks <- dba.peakset(panc_consensus, bRetrieve=TRUE)

# Count reads in the consensus peaks set
panc <- dba.count(panc, peaks=consensus.peaks)

# Normalize the data
panc <- dba.normalize(panc, background=TRUE, normalize=DBA_NORM_NATIVE)

# Perform the differential analysis
panc <- dba.contrast(panc, design="~Condition", reorderMeta=list(Condition="WT"))
panc <- dba.analyze(panc)

# plot a correlation heatmap using the differentially bound sites only
png(file.path(figdir, "heatmap.sample-to-sample.png"), width = 3600, height = 3200, units = "px", pointsize = 6, res=600)
plot(panc, contrast=1)
dev.off()

# retrieve the differentially bound sites
panc.DB <- dba.report(panc)

# profiles on combined conditions
profiles.combined <- dba.plotProfile(panc)

# plot profile heatmap
png(file.path(figdir, 'profile.heatmap.by_condition.png'), width = 6400, height = 6400, units = "px", pointsize = 6, res=600)
dba.plotProfile(profiles.combined)
dev.off()

# Annotate differentially bound sites
# load the knownGene track
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# change the style of the chromosome names from UCSC to Ensembl
seqlevelsStyle(txdb) <- "NCBI"

# prepare the TSS regions, which are defined as the flanking sequence of the TSS sites.
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

# load the DE differential bound sites
de.peaks <- as.data.frame(panc.DB)

# convert to GRanges object
de.peaks.gr <- makeGRangesFromDataFrame(de.peaks, keep.extra.columns=TRUE)

# peak annotation
peakAnno <- annotatePeak(de.peaks.gr, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")

# write to file
write.table(as.data.frame(peakAnno), file=file.path(infodir, 'DE.peaks.annotated.txt'), quote=F, sep='\t', row.names=FALSE)

sessionInfo()

