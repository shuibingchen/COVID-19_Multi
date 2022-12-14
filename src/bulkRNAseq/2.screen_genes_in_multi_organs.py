# 2.screen_genes_in_multi_organs.py
# screen candidate genes that are commonly up-regulated in SARS-CoV-2 infected conditions
# 

import pandas as pd
import os.path
from re import search

workdir = '.'
sourcedir = os.path.join(workdir, 'source')
figdir = os.path.join(workdir, "results", "bulkRNAseq", "figure")
infodir = os.path.join(workdir, "results", "bulkRNAseq", "info")

# collect data
mytable = None
# 1. Airway 48h Infected vs. Mock
for cdt in ['Airway.48h.MOI0.01','Airway.48h.MOI0.1','Airway.48h.MOI1']:
    # read in results
    data = pd.read_table(os.path.join(infodir, 'DE.{}.vs.Airway.48h.Mock.LFC.tsv'.format(cdt)), sep='\t', header=0, usecols=[0,1,2,5,7,8], low_memory=False)
    # rename columns
    data.rename(columns={'baseMean':'{}_baseMean'.format(cdt),
                            'log2FoldChange':'{}_log2FC'.format(cdt),
                            'padj':'{}_padj'.format(cdt)}, inplace=True)
    # append to table
    if mytable == None:
        mytable = data
    else:
        mytable = mytable.merge(data, on=['gene','gene_name','gene_biotype'])
# 2. Alveolar 48h Infected vs. Mock
for cdt in ['Alveolar.48h.MOI0.01','Alveolar.48h.MOI0.1','Alveolar.48h.MOI1']:
    # read in results
    data = pd.read_table(os.path.join(infodir, 'DE.{}.vs.Alveolar.48h.Mock.LFC.tsv'.format(cdt)), sep='\t', header=0, usecols=[0,1,2,5,7,8], low_memory=False)
    # rename columns
    data.rename(columns={'baseMean':'{}_baseMean'.format(cdt),
                            'log2FoldChange':'{}_log2FC'.format(cdt),
                            'padj':'{}_padj'.format(cdt)}, inplace=True)
    # append to table
    if mytable == None:
        mytable = data
    else:
        mytable = mytable.merge(data, on=['gene','gene_name','gene_biotype'])
# 3. Cardiomyocyte 48h Infected vs. Mock
for cdt in ['Cardiomyocyte.48h.MOI0.01','Cardiomyocyte.48h.MOI0.1','Cardiomyocyte.48h.MOI1']:
    # read in results
    data = pd.read_table(os.path.join(infodir, 'DE.{}.vs.Cardiomyocyte.48h.Mock.LFC.tsv'.format(cdt)), sep='\t', header=0, usecols=[0,1,2,5,7,8], low_memory=False)
    # rename columns
    data.rename(columns={'baseMean':'{}_baseMean'.format(cdt),
                            'log2FoldChange':'{}_log2FC'.format(cdt),
                            'padj':'{}_padj'.format(cdt)}, inplace=True)
    # append to table
    if mytable == None:
        mytable = data
    else:
        mytable = mytable.merge(data, on=['gene','gene_name','gene_biotype'])

# save table
mytable.to_csv(os.path.join(infodir, 'DE.48h.summary.txt'), index=False, na_rep='NA', sep='\t')

# screen candidate genes
# cutoff
cut_padj = 0.05
cut_fc = 0.75
cut_bm = 10

# columns to check
cdts = ['Airway.48h.MOI0.01','Airway.48h.MOI0.1','Airway.48h.MOI1',
    'Alveolar.48h.MOI0.01','Alveolar.48h.MOI0.1','Alveolar.48h.MOI1',
    'Cardiomyocyte.48h.MOI0.01','Cardiomyocyte.48h.MOI0.1','Cardiomyocyte.48h.MOI1']

cols_padj = ['{}_padj'.format(x) for x in cdts]
cols_fc = ['{}_log2FC'.format(x) for x in cdts]
cols_bm = ['{}_baseMean'.format(x) for x in cdts]

# perform filtering on genes
# count number of comparisons passing filters
def myFiltCounter(tx, cols_padj, cols_fc, cols_bm, cut_padj, cut_fc, cut_bm):
    nHits = 0
    for padj,fc,bm in zip(cols_padj, cols_fc, cols_bm):
        if tx[padj] < cut_padj and tx[fc] > cut_fc and tx[bm] > cut_bm and tx['gene_biotype'] == 'protein_coding':
            nHits += 1
    return nHits

mytable['Filt'] = mytable.apply(myFiltCounter, axis=1, args=(cols_padj, cols_fc, cols_bm, cut_padj, cut_fc, cut_bm))

# print out candidate genes
print(mytable[mytable['Filt']>=7][['gene_name','Filt']])

print('All Complete!')

