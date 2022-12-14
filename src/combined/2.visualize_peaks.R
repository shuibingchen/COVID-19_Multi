# 2.visualize_peaks.R
# visualize peaks associated with the NR4A1 gene in WT and CIART-/- hPSC-CMs in ATACseq and CUT&RUN assays
# 

library(karyoploteR)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(svglite)

workdir <- '.'
sourcedir <- file.path(workdir, "data")

# bigwig files for CUT&RUN samples can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202966
cutrun.bwdir <- file.path(sourcedir, "cutrun", "bigwig")

# bigwig files for ATACseq samples can be downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202965
atac.bwdir <- file.path(sourcedir, "atac", "bigwig")

figdir <- file.path(workdir, "results", "combined", "figure")
infodir <- file.path(workdir, "results", "combined", "info")

cutrun.bw.files <- file.path(cutrun.bwdir, c('WT_1.bw','WT_2.bw','IgG_Control.bw'))
atac.bw.files <- file.path(atac.bwdir, c('CM_WT_1.bw', 'CM_WT_2.bw', 'CM_WT_3.bw', 'CM_CIART_1.bw', 'CM_CIART_2.bw', 'CM_CIART_3.bw'))

# function to plot peaks
my.plot.peaks <- function(ref.genome="hg19", bsgenome=BSgenome.Hsapiens.UCSC.hg19, knowngene=TxDb.Hsapiens.UCSC.hg19.knownGene, 
    seq.level.style='NCBI', gene="NR4A1", flank=2000, ymax="visible.region", color="#0000FF", big.wig.files=c(atac.bw.files)){
  # Plot parameters, only to look better
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  pp$data1outmargin <- 0

  # Get transcrupts annotation to get gene regions
  tssAnnot <- ELMER::getTSS(genome = ref.genome)
  seqlevelsStyle(tssAnnot) <- seq.level.style
  tssAnnot <- tssAnnot[tssAnnot$external_gene_name == gene]
  
  bsgenome <- bsgenome
  seqlevelsStyle(bsgenome) <- seq.level.style
  
  knowngene <- knowngene
  seqlevelsStyle(knowngene) <- seq.level.style
  
  # plot will be at the gene range +- 2 Kb
  gene.region <- range(c(tssAnnot)) + flank
  
  # Start by plotting gene tracks
  # chromosome
  kp <- plotKaryotype(zoom = gene.region,genome = bsgenome, cex = 0.5, plot.params = pp)
  genes.data <- makeGenesDataFromTxDb(knowngene,
                                      karyoplot = kp,
                                      plot.transcripts = TRUE, 
                                      plot.transcripts.structure = TRUE)
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  
  # gene
  kp <- plotKaryotype(zoom = gene.region,genome = bsgenome, cex = 0.5, plot.params = pp)
  kpAddBaseNumbers(kp, tick.dist = 20000, minor.tick.dist = 5000,
                   add.units = TRUE, cex = 0.4, tick.len = 3)
  kpPlotGenes(kp, data = genes.data, r0 = 0, r1 = 0.25, gene.name.cex = 0.5)
  
  # reserve area to plot the bigwig files
  big.wig.files <- big.wig.files
  out.at <- autotrack(1:length(big.wig.files), 
                      length(big.wig.files), 
                      margin = 0.3, 
                      r0 = 0.3,
                      r1 = 1)
  
  for(i in seq_len(length(big.wig.files))) {
    bigwig.file <- big.wig.files[i]
    
    # Define where the track will be plotted
    # autotrack will simple get the reserved space (from out.at$r0 up to out.at$r1)
    # and split in equal sizes for each bigwifile, i the index, will control which 
    # one is being plotted
    at <- autotrack(i, length(big.wig.files), r0 = out.at$r0, r1 = out.at$r1, margin = 0.2)
    
    # Plot bigwig
    kp <- kpPlotBigWig(kp, 
                       data = bigwig.file, 
                       ymax = ymax,
                       r0 = at$r0, 
                       #col = ifelse(grepl("ATACseq",bigwig.file),"#0000FF","#FF0000"),
                       col = color,
                       r1 = at$r1)
    computed.ymax <- ceiling(kp$latest.plot$computed.values$ymax)
    
    # Add track axis
    kpAxis(kp, 
           ymin = 0, 
           ymax = computed.ymax, 
           numticks = 2,
           r0 = at$r0, 
           r1 = at$r1,
           cex = 0.5)
    
    # Add track label
    kpAddLabels(kp, 
                ##labels = ifelse(grepl("ATACseq",bigwig.file),"ATACseq","CUT&RUN"),
                labels = gsub(".BPM.bw", "", basename(bigwig.file)),
                r0 = at$r0, 
                r1 = at$r1, 
                cex = 0.5, 
                label.margin = 0.01)
  }
}

# plot ATACseq peaks
svg(file.path(figdir, 'ATACseq.peaks.NR4A1.svg'), width=5, height=5)
my.plot.peaks(ref.genome="hg19", bsgenome=BSgenome.Hsapiens.UCSC.hg19, knowngene=TxDb.Hsapiens.UCSC.hg19.knownGene, seq.level.style='NCBI', gene="NR4A1", flank=10000, ymax=3, color="#386cb0", big.wig.files=c(atac.bw.files))
dev.off()

# plot CUT&RUN peaks
svg(file.path(figdir, 'CUTRUN.peaks.NR4A1.svg'), width=5, height=3.5)
my.plot.peaks(ref.genome="hg19", bsgenome=BSgenome.Hsapiens.UCSC.hg19, knowngene=TxDb.Hsapiens.UCSC.hg19.knownGene, seq.level.style='UCSC', gene="NR4A1", flank=10000, ymax=4, color="#d95f02", big.wig.files=c(cutrun.bw.files))
dev.off()

sessionInfo()

