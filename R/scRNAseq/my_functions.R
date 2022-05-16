# my_functions.R
# customized functions for processing data and plotting with Seurat3
# Author: Tuo Zhang
# Date: 10/09/2019
# Version: 1.0
# NEW: first version
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
library(scater, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)
library(cowplot, quietly=T, warn.conflicts=F)
library(R.utils, quietly=T, warn.conflicts=F)
library(SingleCellExperiment, quietly=T, warn.conflicts=F)
library(tidyr, quietly=T, warn.conflicts=F)
library(tibble, quietly=T, warn.conflicts=F)
library(ggrepel, quietly=T, warn.conflicts=F)
#library(EnsDb.Mmusculus.v79, quietly=T)
#library(EnsDb.Hsapiens.v86, quietly=T)

# read in raw counts data from a given sample
my.Read10X <- function(tcountdir, tprefix=NULL){
  # directory exists?
  if (! dir.exists(tcountdir)){
    print(paste("input raw counts folder does NOT exist:", tcountdir, sep=" "))
    return(NULL)
  }
  # file exists?
  tmat.file <- paste(tcountdir, "matrix.mtx.gz", sep="/")
  tfnames.file <- paste(tcountdir, "features.tsv.gz", sep="/")
  tbnames.file <- paste(tcountdir, "barcodes.tsv.gz", sep="/")
  for (tf in c(tmat.file, tfnames.file, tbnames.file)){
    if (! file.exists(tf)){
      print(paste("input file does NOT exist:", tf))
      return(NULL)
    }
  }
  # extract counts matrix
  cat(paste("Loading UMI counts table from", tcountdir, "..."))
  tmat <- readMM(paste(tcountdir, "matrix.mtx.gz", sep="/"))
  tfnames <- read.delim(paste(tcountdir, "features.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  tbnames <- read.delim(paste(tcountdir, "barcodes.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  # update column names (cell ids)
  if (is.null(tprefix)){
    colnames(tmat) = tbnames$V1
  }
  else{
    colnames(tmat) = paste(tprefix, tbnames$V1, sep="_")
  }
  #rownames(tmat) = tfnames$V1
  # replace rowname (Ensembl id) by gene symbol
  # in case gene symbol is not unique, append the _EnsemblID after it
  # missing gene symbol will be replaced by EnsemblID
  #tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="GENENAME", keytype="GENEID")
  ##tsymbols <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(tmat), column="SYMBOL", keytype="GENEID")
  rownames(tmat) <- uniquifyFeatureNames(ID=tfnames$V1, names=tfnames$V2)
  cat(" done.","\n")
  return(tmat)
}

# merge raw read counts table collected from multiple samples
my.MergeMatrix <- function(tmats){
  cat("Merge raw UMI counts ")
  tfunc <- function(x,y){
    if (class(x) != "data.frame")
      x <- as.data.frame(as.matrix(x))
    if (class(y) != "data.frame")
      y <- as.data.frame(as.matrix(y))
    tres <- merge(x, y, by=0, all=T)
    rownames(tres) <- tres$Row.names
    tres <- tres[,-1]
    cat(".")
    return(tres)
  }
  tmerged <- Reduce(f=tfunc, x=tmats)
  # fill na with 0
  tmerged[is.na(tmerged)] <- 0
  tmerged <- as(as.matrix(tmerged), "dgCMatrix")
  cat(" done.")
  return(tmerged)
}

# plot cells colored by an annotation (DimPlot)
myDimPlot <- function(tobj, treduct, tcate, torder=NULL, tsuffix, tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit=FALSE, txlim=NULL, tylim=NULL,
                      tncol=4, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      missed.categroy <- setdiff(unique(tdataToPlot$Category), torder)
      if (length(missed.categroy) > 0){
        tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=missed.categroy))
      }
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, color=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      missed.categroy <- setdiff(unique(tdataToPlot$Category), torder)
      if (length(missed.categroy) > 0){
        tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=missed.categroy))
      }
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, color=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

# highlight expression of a set of genes (FeaturePlot)
MyFeaturePlot <- function(tobj, tgenes, tcells=NULL, tassay="RNA", treduction.name="umap", tsort=FALSE, txlim=NULL, tylim=NULL, tbreaks=NULL, tlimits=NULL, tlowcolor='gray80', thighcolor='red2', tncol=2, tlegend=NULL, tptsize=2, talpha=0.7, twidth=15, theight=12.5, tunits="in", tres=300){
  # genes valid?
  tgenes.valid <- intersect(tgenes, rownames(tobj))
  if (is.null(tgenes.valid)){
    cat("No valid genes found, do nothing!")
    return(NULL)
  }
  # assay valid?
  if (! tassay %in% names(tobj)){
    cat(paste("Not a valid assay:",tassay,sep=" "))
    return(NULL)
  }
  # extract gene expression
  texp <- as.matrix(tobj[[tassay]]@data[tgenes.valid, ,drop=F])
  # get coordinates
  tvars <- c("UMAP_1","UMAP_2")
  if (treduction.name == "tsne"){
    tvars <- c("tSNE_1","tSNE_2")
  }
  tdata <- FetchData(object=tobj, vars=tvars)
  colnames(tdata) <- c("X","Y")
  # plot
  tplots <- list()
  tk <- 1
  for (tgene in tgenes){
    # merge data for plotting
    tdata.merged <- merge(tdata, t(texp[tgene,,drop=F]), by=0, all=T)
    rownames(tdata.merged) <- tdata.merged$Row.names
    tdata.merged <- tdata.merged[,-1]
    colnames(tdata.merged) <- c("X","Y","Expression")
    # subset cells?
    if (! is.null(tcells)){
      tdata.merged <- tdata.merged[tcells,,drop=F]
    }
    # reorder cells by expression of the given gene
    if (tsort){
      tdata.merged <- tdata.merged[with(tdata.merged, order(Expression)),]
    }
    # plot
    if (max(tdata.merged$Expression) > 0){ # expressed in at least one cell
      # plot (rename x and y axis)
      tg <- ggplot(tdata.merged, aes(x=X, y=Y, color=Expression))
      tg <- tg + geom_point(shape=19, size=tptsize, alpha=talpha)
      if (! is.null(tbreaks)){
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor, breaks=tbreaks, limits=tlimits)
      } else {
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor)
      }
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    } else { # no expressions at all
      tg <- ggplot(tdata.merged, aes(x=X, y=Y))
      tg <- tg + geom_point(color="gray80", shape=19, size=tptsize, alpha=talpha)
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    }
  }
  # combine plots with Seurat::CombinePlots
  tcombined <- CombinePlots(tplots, ncol=tncol, legend=tlegend)
  return(tcombined)
  # ggsave(paste(toutdir, paste("gene","exp",tsuffix,toupper(treduction.name),"png",sep="."),sep="/"), height=theight, width=twidth, units=tunits, dpi=tres)
  # ggsave(paste(toutdir, paste("gene","exp",tsuffix,toupper(treduction.name),"svg",sep="."),sep="/"), height=theight, width=twidth, units=tunits)
}

# running original DE test, without considering conservations in each sample
myDETest <- function(tobj, tk, tmethod, tplot, tinfodir, tfigdir, tassay=NULL, tsuf=NULL){
  # DE test
  tmarkers <- FindMarkers(object=tobj, ident.1=tk, min.pct=0.25, test.use=tmethod, assay=tassay)
  # file name
  tdesp <- ".origin.markers"
  if (! is.null(tsuf)){
    tdesp <- paste("",tsuf,"origin","markers",sep=".")
  }
  # write marker genes to file
  write.table(as.data.frame(tmarkers), file=paste(tinfodir, paste("C", tk, ".", tmethod, tdesp, ".txt", sep=""), sep="/"), quote=FALSE, na="", sep="\t", col.names=NA)
  # select top 4 positive markers
  genes.viz <- head(rownames(subset(tmarkers, avg_logFC>0)),4)
  print(genes.viz)
  # visualize markers with a violin plot
  if(tplot){
    print("violin plot 1")
    png(file=paste(tfigdir, paste("VlnPlot.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""), sep="/"), width=6400, height=6400, units="px", pointsize=6, res=600)
    print(VlnPlot(object=tobj, features=genes.viz, ncol=2))
    dev.off()
    # view which cells express a given gene (red is high expression) on a tSNE plot
    # print("tSNE plot")
    # png(file=paste(tfigdir, paste("gene.exp.tSNE.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""),sep="/"), width=4000, height=3200, units="px", pointsize=6, res=600)
    # print(FeaturePlot(object=tobj, features=genes.viz, reduction="tsne", pt.size=1, cols=c("grey","red")))
    # dev.off()
    # # view which cells express a given gene (red is high expression) on a UMAP plot
    print("UMAP plot")
    png(file=paste(tfigdir, paste("gene.exp.UMAP.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""),sep="/"), width=4000, height=3200, units="px", pointsize=6, res=600)
    print(FeaturePlot(object=tobj, features=genes.viz, reduction="umap", pt.size=1, cols=c("grey","red")))
    dev.off()
  }
}

# stacked bar plot showing the percentage of cells per cluster in each group
myStackBarCellComposition <- function(tobj, tcluster_by='ident', tcluster_order, tgroup_by='Name', tgroup_order=NULL, tcells=NULL, tanncolor=NULL, ttlsize=28, ttxsize=24, tltsize=18){
  # calculate percentage of cells per cluster in each group
  tpercent <- FetchData(tobj, vars=c(tgroup_by, tcluster_by), cells=tcells) %>% dplyr::rename('Cluster'=tcluster_by, 'Group'=tgroup_by) %>%
    rownames_to_column('CellID') %>% group_by(Group, Cluster) %>% summarise_at('CellID', length) %>%
    tidyr::complete(Cluster=factor(tcluster_order), fill=list(CellID=0)) %>%
    group_by(Group) %>% mutate_at('CellID', function(x) { round(x/sum(x)*100,2) })
  colnames(tpercent) <- c("Group", "Cluster", "Percent")
  # reorder group/cluster if provided
  if (! is.null(tcluster_order)){
    tpercent$Cluster <- factor(tpercent$Cluster, levels=tcluster_order)
  }
  if (! is.null(tgroup_order)){
    tpercent$Group <- factor(tpercent$Group, levels=tgroup_order)
  }
  # plot
  tg <- ggplot(tpercent, aes(x=Group, y=Percent, fill=Cluster))
  tg <- tg + geom_bar(stat='identity')
  if (! is.null(tanncolor)){
    tg <- tg + scale_fill_manual(values=tanncolor)
  }
  tg <- tg + ylab('Percentage of cells')
  tg <- tg + theme_bw()
  tg <- tg + theme(axis.title.x=element_blank(), axis.title.y=element_text(color='black', size=ttlsize))
  tg <- tg + theme(axis.text.x=element_text(color='black', size=ttxsize, angle=60, hjust=1), axis.text.y=element_text(color='black', size=ttxsize))
  tg <- tg + theme(legend.text=element_text(size=tltsize), legend.title=element_blank())
  return(tg)
}

# calculate the fraction of overlapping marker genes in the reference data
calOverlapFrac <- function(tref, tgenes){
  return(length(intersect(tref, tgenes)) / length(tref))
}

# draw heatmap showing the fraction of overlapping marker genes between reference and our data
plotFracOverlap <- function(trefList, tqryList, tinfodir, tfigdir, toutprefix, tsvg=FALSE){
  # calculate fraction matrix
  fracMat <- matrix(nrow=length(tqryList), ncol=length(trefList))
  colnames(fracMat) <- names(trefList)
  rownames(fracMat) <- names(tqryList)
  for (i in names(tqryList)){
    for (j in names(trefList)){
      fracMat[i,j] <- calOverlapFrac(trefList[[j]], tqryList[[i]])
    }
  }
  # save fraction matrix
  write.table(fracMat, file.path(tinfodir, paste0(toutprefix, ".txt")), quote=F, sep='\t', col.names=NA)
  # plot overlapping fraction matrix
  # normalized fraction
  pheatmap(fracMat, color=colorRampPalette(c("#4575B4", "#f7f7f7", "#D73027"))(100), cluster_rows=F, cluster_cols=F, 
           scale='row', show_rownames=T, show_colnames=T, fontsize=10,
           filename=file.path(tfigdir, paste0(toutprefix, '.rescaled.png')), width=7, height=3)
  if (tsvg){
    svg(file.path(tfigdir, paste0(toutprefix, '.rescaled.svg')), width=7, height=3)
    pheatmap(fracMat, color=colorRampPalette(c("#4575B4", "#f7f7f7", "#D73027"))(100), cluster_rows=F, cluster_cols=F, 
             scale='row', show_rownames=T, show_colnames=T, fontsize=10)
    dev.off()
  }
}

# score each cell based on pre-selected ref marker genes, and plot heatmap/UMAP
scoreByRefMarkers <- function(tobj, trefList, trefPattern='^Ref_', tinfodir, tfigdir, toutprefix, tsvg=FALSE){
  # run AddModuleScore to score each reference marker set
  tenrich.name <- "REF"
  tp <- AddModuleScore(tobj, features=trefList, pool=NULL, nbin=24, ctrl=min(vapply(X=trefList, FUN=length, FUN.VALUE=numeric(length=1))), k=FALSE, assay=NULL, name=tenrich.name, seed=1)
  # collect reference set score per cell
  tref.scores <- FetchData(tp, vars=c('ident', grep(tenrich.name, colnames(tp@meta.data), value=T)))
  colnames(tref.scores) <- c('ident', names(trefList))
  # write per-cell reference set scores to file
  write.table(tref.scores, file=file.path(tinfodir,paste(toutprefix,"txt",sep=".")), quote=F, sep="\t", row.names=T, col.names=NA)
  
  # assign each cell a reference cluster if possible
  assignClust <- function(tx){
    if (all(tx < 0)){
      return('Unknown')
    } else if (length(which(tx == max(tx))) > 1){
      return('Draw')
    } else {
      return(names(tx)[which(tx == max(tx))])
    }
  }
  ref.assign <- apply(tref.scores[,grep(trefPattern, colnames(tref.scores), value=T)], MARGIN=1, FUN=assignClust)
  #table(ref.assign)
  # normalize the scores per cell
  tref.scores.norm <- as.data.frame(t(scale(t(tref.scores[,grep(trefPattern, colnames(tref.scores), value=T)]))))
  colnames(tref.scores.norm) <- paste(tenrich.name, 1:length(trefList),'Norm',sep='_')
  # add to seurat object
  tp[['REF.assign']] <- ref.assign
  tp$REF.assign <- factor(tp$REF.assign, levels=c(names(trefList), 'Unknown'))
  tp[[colnames(tref.scores.norm)]] <- tref.scores.norm
  
  # UMAP showing normalized scores
  tdata <- FetchData(tp, vars=c('UMAP_1','UMAP_2', paste(tenrich.name, 1:length(trefList),'Norm',sep='_')))
  
  for (k in 1:length(trefList)){
    tclust <- paste(tenrich.name,k,'Norm',sep='_')
    print(tclust)
    g <- ggplot(tdata, aes_string(x='UMAP_1', y='UMAP_2', color=tclust))
    g <- g + geom_point(shape=19, size=1.5, alpha=0.7)
    g <- g + scale_color_gradient2(low="#4575B4", mid="#f7f7f7", high="#D73027")
    g <- g + ggtitle(tclust)
    g <- g + xlab("UMAP 1") + ylab("UMAP 2")
    g <- g + theme_bw()
    g <- g + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
    g <- g + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
    ggsave(file.path(tfigdir, paste('UMAPPlot',toutprefix, tclust,'png',sep='.')), plot=g, width=10, height=8, dpi=300)
    if (tsvg){
      ggsave(file.path(tfigdir, paste('UMAPPlot',toutprefix, tclust,'svg',sep='.')), plot=g, width=8, height=6.4)
    }
  }
} 
