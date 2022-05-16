from statsmodels.stats import multitest
from scipy import stats
import pandas as pd
import numpy as np

mop_df = pd.read_csv('../outputs/gene_cell/Merfish_MOp_peripheral.csv')
mop_df['sample_id'] = 'm'+mop_df['mouse'].astype(str)+'s'+mop_df['sample'].astype(str)

#filter SRRS results to
#1. Drop gene/cells which have fewer than 5 gene spots of interest
#2. Drop gene/cell-type/samples with fewer than 20 cells
mop_df = mop_df[mop_df['num_gene_spots'].ge(5)]
mop_df = mop_df.groupby(['sample_id','gene','annotation']).filter(lambda g: len(g) >= 20)

#Calculate z from Lyapunov CLT for each gene in each sample
gb_cols = ['sample_id','gene','annotation']

mop_agg_df = mop_df.groupby(gb_cols).agg(
    num_cells = ('cell_id','nunique'),
    med_gene_spots = ('num_gene_spots','median'),
    med_spots = ('num_spots','median'),
    med_score = ('score','median'),
    score_sum = ('score','sum'),
    var_sum = ('variance','sum'),
).reset_index()

mop_agg_df['z'] = mop_agg_df['score_sum']/np.sqrt(mop_agg_df['var_sum'])

#Calculate two-sided p and BH correct ps
p_onesided = stats.norm.cdf(mop_agg_df['z'])
mop_agg_df['p'] = 2*np.minimum(p_onesided, 1-p_onesided)

alpha = 0.05

_,adj_p,_,_ = multitest.multipletests(
    mop_agg_df['p'],
    alpha = alpha,
    method = 'fdr_bh',
)
mop_agg_df['bh_p'] = adj_p

mop_agg_df.to_csv('../outputs/gene_ontology/Merfish_MOp_peripheral.csv',index=False)

