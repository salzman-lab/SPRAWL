from statsmodels.stats import multitest
from scipy import stats
import pandas as pd
import numpy as np

viz_df = pd.read_csv('../outputs/gene_cell/Viz_central.csv')
viz_df = viz_df.rename(columns={'annotation':'ontology'})
#viz_df['sample_id'] = 's'+viz_df['sample'].astype(str)+'r'+viz_df['replicate'].astype(str)

#filter SRRS results to
#1. Drop gene/cells which have fewer than 5 gene spots of interest
#2. Drop gene/cell-type/samples with fewer than 20 cells
viz_df = viz_df[viz_df['num_gene_spots'].ge(5)]
viz_df = viz_df.groupby(['sample_id','gene','ontology']).filter(lambda g: len(g) >= 20)

#Calculate z from Lyapunov CLT for each gene in each sample
gb_cols = ['sample_id','gene','ontology']

viz_agg_df = viz_df.groupby(gb_cols).agg(
    num_cells = ('cell_id','nunique'),
    med_gene_spots = ('num_gene_spots','median'),
    med_spots = ('num_spots','median'),
    med_score = ('score','median'),
    score_sum = ('score','sum'),
    var_sum = ('variance','sum'),
).reset_index()

viz_agg_df['z'] = viz_agg_df['score_sum']/np.sqrt(viz_agg_df['var_sum'])

#Calculate two-sided p and BH correct ps
p_onesided = stats.norm.cdf(viz_agg_df['z'])
viz_agg_df['p'] = 2*np.minimum(p_onesided, 1-p_onesided)

alpha = 0.05

_,adj_p,_,_ = multitest.multipletests(
    viz_agg_df['p'],
    alpha = alpha,
    method = 'fdr_bh',
)
viz_agg_df['bh_p'] = adj_p

viz_agg_df.to_csv('../outputs/gene_ontology/Viz_central.csv',index=False)


