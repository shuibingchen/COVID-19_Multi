# 3.plot_candidate_genes.R
# Generate 3D plots and highlight candidate multi-organs SARS-CoV-2 infections RNAseq samples
# 

library(tidyverse)
library(plotly)
library(htmlwidgets)
library(circlize)

workdir <- '.'
figdir <- file.path(workdir, "results", "bulkRNAseq", "figure")
infodir <- file.path(workdir, "results", "bulkRNAseq", "info")

# read in DE table
data <- read.table(file.path(infodir, 'DE.48h.summary.txt'), header=T, check.names=F, stringsAsFactors=F, sep='\t')

# pre-filter on the table:
# 1. protein coding genes only
# 2. replace NAs (padj columns) with 1 (insignificance)
data <- data %>% filter(gene_biotype == 'protein_coding') %>% 
  replace(is.na(.), 1)

# pre-label the candidate genes
candidates <- data %>% filter(Filt >= 7) %>% select(c('gene','gene_name'))
data$label <- with(data, ifelse(gene %in% candidates$gene, gene_name, ''))

# make 3D scatter plot per organ per MOI
# In each plot, one point indicates one gene (Ensembl id); three dimensions represent the three metrics used for screening:
# - logBaseMean: log10(base Mean + 1)
# - logFoldChange
# - logPadj: -log10(adjusted p-value)
# 
# Highlight points (genes) passing filter:
# - basemean > 10
# - logFoldChange > 0.75
# - padj < 0.05
# 
# for cases with padj == 0, replace 0 by 1e-320
# 

# function to make plot for a given cohort (no label/color candidate genes)
my.plot.metrics.no.label <- function(data, cohort, bm.cut, lfc.cut, padj.cut, colors=c('gray65','#cab2d6')){
  # related columns
  my.columns <- c('gene', 'gene_name', 'label', paste(cohort, 'baseMean', sep='_'),
                  paste(cohort, 'log2FC', sep='_'),
                  paste(cohort, 'padj', sep='_'))
  # rename columns
  my.rename <- c(paste(cohort, 'baseMean', sep='_'), 
                 paste(cohort, 'log2FC', sep='_'),
                 paste(cohort, 'padj', sep='_'))
  names(my.rename) <- c('baseMean', 'log2FC', 'padj')
  # clean data for plotting
  data.to.plot <- data %>% select(my.columns) %>% 
    rename(my.rename) %>% 
    mutate(padj.fix=ifelse(padj==0, 1e-320, padj)) %>% 
    mutate(logBaseMean=log10(baseMean+1), 
           logPadj=-log10(padj.fix),
           candidate=ifelse(baseMean > get("bm.cut", envir=.env) & 
                              log2FC > get("lfc.cut", envir=.env) & 
                              padj < get("padj.cut", envir=.env), 
                            'passed','filtered')) %>% 
    mutate(candidate=factor(candidate, levels=c('filtered','passed'))) %>% 
    arrange(desc(candidate))
  # variables for x,y,z
  x <- 'logBaseMean'
  y <- 'log2FC'
  z <- 'logPadj'
  # make plot
  fig <- plot_ly(data.to.plot, x=as.formula(paste0('~',x)),
                 y=as.formula(paste0('~',y)),
                 z=as.formula(paste0('~',z)),
                 color=~candidate,
                 colors=colors) %>%
    add_markers(alpha=0.9) %>% 
    layout(title=paste(cohort, "label", sep='-'),
           scene=list(xaxis=list(title='logBaseMean'),
                      yaxis=list(title='logFoldChange'),
                      zaxis=list(title='logAdjustedPvalue'),
                      aspectmode='cube'))
    ##add_trace(x=3.680417, y=1.847951, z=81.141838, type="scatter3d", text="HIVEP2", mode="text")
  return(fig)
}

# cutoff
bm.cut <- 10
lfc.cut <- 0.75
padj.cut <- 0.05

# plot
for (cohort in gsub('_padj','', grep('padj', colnames(data), value=T))){
  print(my.plot.metrics.no.label(data, cohort, bm.cut, lfc.cut, padj.cut, colors=c('gray65','#984ea3')))
}

sessionInfo()

