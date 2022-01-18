from statsmodels.stats import multitest
import pandas as pd
import numpy as np
from scipy import stats

import glob
import sys
import os

#Read the gene_cell level table
gene_cell_df = pd.read_csv('../outputs/gene_cell/Vizgen_Brainmap_peripheral.csv')

#Calculate z from Lyapunov CLT for each gene in each sample
gb_cols = ['mouse','sample','replicate','gene']

gene_df = gene_cell_df.groupby(gb_cols).agg(
    num_cells = ('cell_id','nunique'),
    num_annotations = ('annotation','nunique'),
    med_gene_spots = ('num_gene_spots','median'),
    med_spots = ('num_spots','median'),
    med_score = ('score','median'),
    score_sum = ('score','sum'),
    var_sum = ('variance','sum'),
).reset_index()

gene_df['z'] = gene_df['score_sum']/np.sqrt(gene_df['var_sum'])

#Calculate two-sided p and BH correct ps
p_onesided = stats.norm.cdf(gene_df['z'])
gene_df['p'] = 2*np.minimum(p_onesided, 1-p_onesided)

alpha = 0.05

_,adj_p,_,_ = multitest.multipletests(
    gene_df['p'],
    alpha = alpha,
    method = 'fdr_bh',
)
gene_df['bh_p'] = adj_p

gene_df.to_csv('../outputs/gene/Vizgen_Brainmap_peripheral.csv',index=False)

