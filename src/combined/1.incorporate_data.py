# 1.incorporate_data.py
# incorporate data from CUT&RUN, ATACseq and RNAseq, to screen peaks(genes) that are bound by CIART and 
# are significantly changed at both transcriptional and chromatin accessibility levels,
# to identify the downstream targets regulated by CIART in hPSC-CMs
# 

import sys
import os.path
import pandas as pd
import numpy as np
import re

workdir = "."
sourcedir = os.path.join(workdir, "source")
figdir = os.path.join(workdir, "results", "combined", "figure")
infodir = os.path.join(workdir, "results", "combined", "info")

cutrundir = os.path.join(workdir, "results", "cutrun", "info")
atacdir = os.path.join(workdir, "results", "ATACseq", "info")
rnaseqdir = os.path.join(workdir, "results", "bulkRNAseq", "info")

cutrun_file = os.path.join(cutrundir, "annotated.peaks.candidate.txt")
atacseq_file = os.path.join(atacdir, "DE.peaks.annotated.txt")
rnaseq_file = os.path.join(rnaseqdir, "DE.Cardiomyocyte.CIART.Mock.vs.Cardiomyocyte.WT.Mock.LFC.tsv")

# Load ATACseq data and summarize it at gene level
# read in ATACseq data
atac = pd.read_table(atacseq_file, header=0, sep='\t', low_memory=False)

# clean data by: 
# 1) removing entries for which the Ensembl ids are unavailable
# 2) condense entries that targets the same gene
def series2str(ts, sep=';'):
    """pd.series to ';' split string"""
    return sep.join(['{}'.format(x) for x in ts.tolist()])

def condense(td):
    """given a groupby data block, condense information to a single entry"""
    # create a series for storing peak coordinates
    coords = td.apply(lambda x:'{}:{}-{}'.format(x['seqnames'],x['start'],x['end']), axis=1)
    # hits promoter regions?
    promoters = td['annotation'].str.contains('Promoter')
    # more open/close on promoter region?
    promoters_open = promoters_close = False
    if promoters.any():
        promoters_open = (td[promoters]['Fold'] > 0).any()
        promoters_close = (td[promoters]['Fold'] < 0).any()
    # condense information
    return pd.Series([series2str(coords), series2str(td['annotation']), 
                      series2str(td['Fold']), series2str(td['distanceToTSS']), 
                      promoters.any(), promoters_open, promoters_close], 
                     index=['ATAC_coordinates','ATAC_annotations','ATAC_fold','ATAC_distToTSS',
                            'ATAC_promoter','ATAC_promoter_open','ATAC_promoter_close'])

atac_compact = atac[~atac['ENSEMBL'].isna()].groupby('ENSEMBL').apply(condense).reset_index()

# Load CUT&RUN data
# read in Cut&Run data
cutrun = pd.read_table(cutrun_file, header=0, sep='\t', low_memory=False)

# simplify the genomic location annotation
cutrun['simple.annotation'] = cutrun['annotation'].apply(lambda x:x.split(' (')[0])

# Load RNAseq data
# read in RNAseq data
rna = pd.read_table(rnaseq_file, sep='\t', header=0, usecols=[0,1,2,5,7,8], low_memory=False)

# clean DE results:
# 1) use significantly DE genes only (padj < 0.1)
# 2) rename needed columns for merging table purpose
rna = rna[rna['padj'] < 0.1].rename(columns={'gene':'ENSEMBL','baseMean':'RNA_baseMean','log2FoldChange':'RNA_logFC','padj':'RNA_padj'})

# Merge data tables
# merge RNAseq to the Cut&Run base table
rna_cols = ['ENSEMBL','RNA_baseMean','RNA_logFC','RNA_padj']
mytable = pd.merge(cutrun, rna[rna_cols], how='left', on='ENSEMBL')

# fill NaNs with -1, 1 or 0
mytable['RNA_baseMean'] = mytable['RNA_baseMean'].fillna(value=-1)
mytable['RNA_logFC'] = mytable['RNA_logFC'].fillna(value=0)
mytable['RNA_padj'] = mytable['RNA_padj'].fillna(value=1)

# merge ATACseq to the merged table
mytable = pd.merge(mytable, atac_compact, how='left', on='ENSEMBL')

# fill NaNs with False
mytable['ATAC_promoter'] = mytable['ATAC_promoter'].fillna(value=False)
mytable['ATAC_promoter_open'] = mytable['ATAC_promoter_open'].fillna(value=False)
mytable['ATAC_promoter_close'] = mytable['ATAC_promoter_close'].fillna(value=False)

# perform filtering on CUT&RUN peaks based on the other two assays
# 1) ATACseq: changes in chromatin accessibility of the promoter region of the corresponding gene
# 2) RNAseq: expression changes significantly on the corresponding gene (padj < 0.1 + |logFC| > 1)
filt_1 = mytable['ATAC_promoter']
filt_2 = (mytable['RNA_padj'] < 0.1) & ((mytable['RNA_logFC'] > 1) | (mytable['RNA_logFC'] < -1))

# number of peaks before and after filtering steps
print(mytable.shape[0])
print((filt_1).sum())
print((filt_1 & filt_2).sum())

# how many genes hit by these peaks?
len(mytable[filt_1 & filt_2]['SYMBOL'].unique())

# write filltered peaks to file
filtered_table_file = os.path.join(infodir, 'filtered_peaks.xlsx')
mytable[filt_1 & filt_2].to_excel(filtered_table_file, index=False, engine='xlsxwriter')

print('All Complete!')

