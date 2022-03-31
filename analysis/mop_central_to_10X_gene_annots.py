from statsmodels.stats import multitest
from scipy import stats
import pandas as pd
import numpy as np

mop_df = pd.read_csv('../outputs/gene_cell/MOp_central.csv')
mop_df = mop_df.rename(columns={'annotation':'ontology'})
mop_df['sample_id'] = 'm'+mop_df['mouse'].astype(str)+'s'+mop_df['sample'].astype(str)

#filter SRRS results to
#1. Drop gene/cells which have fewer than 5 gene spots of interest
#2. Drop gene/cell-type/samples with fewer than 20 cells
mop_df = mop_df[mop_df['num_gene_spots'].ge(5)]
mop_df = mop_df.groupby(['sample_id','gene','ontology']).filter(lambda g: len(g) >= 20)

#convert MOp annotations to match ReadZs annotations

mop_to_10x_ann_map = {
    'Astro_1':'Astro','Astro_2':'Astro','Astro_3':'Astro',
    'Endo':'Endo',
    'L23_IT_1':'L2/3 IT','L23_IT_2':'L2/3 IT','L23_IT_3':'L2/3 IT','L23_IT_4':'L2/3 IT','L23_IT_5':'L2/3 IT',
    #'L45_IT_1','L45_IT_2','L45_IT_3','L45_IT_4','L45_IT_5','L45_IT_SSp_1','L45_IT_SSp_2',
    'L56_NP_1':'L5/6 NP','L56_NP_2':'L5/6 NP',
    'L5_ET_1':'L5 ET','L5_ET_2':'L5 ET','L5_ET_3':'L5 ET','L5_ET_4':'L5 ET','L5_ET_5':'L5 ET','L5_IT_1':'L5 IT','L5_IT_2':'L5 IT','L5_IT_3':'L5 IT','L5_IT_4':'L5 IT',
    'L6_CT_1':'L6 CT','L6_CT_2':'L6 CT','L6_CT_3':'L6 CT','L6_CT_4':'L6 CT','L6_CT_5':'L6 CT','L6_CT_6':'L6 CT','L6_CT_7':'L6 CT','L6_CT_8':'L6 CT','L6_CT_9':'L6 CT',
    'L6_IT_1':'L6 IT','L6_IT_2':'L6 IT','L6_IT_3':'L6 IT',
    'L6_IT_Car3':'L6 IT Car3',
    'L6b_1':'L6b','L6b_2':'L6b','L6b_3':'L6b','Lamp5_1':'Lamp5',
    'Lamp5_2':'Lamp5','Lamp5_3':'Lamp5','Lamp5_4':'Lamp5','Lamp5_5':'Lamp5','Lamp5_6':'Lamp5','Lamp5_7':'Lamp5','Lamp5_8':'Lamp5','Lamp5_9':'Lamp5',
    #'Micro_1','Micro_2',
    'OPC':'OPC',
    'Oligo_1':'Oligo','Oligo_2':'Oligo','Oligo_3':'Oligo',
    #'PVM',
    #'Peri',
    'Pvalb_1':'Pvalb','Pvalb_10':'Pvalb','Pvalb_11':'Pvalb','Pvalb_12':'Pvalb','Pvalb_2':'Pvalb','Pvalb_3':'Pvalb','Pvalb_4':'Pvalb','Pvalb_5':'Pvalb','Pvalb_6':'Pvalb','Pvalb_7':'Pvalb','Pvalb_8':'Pvalb','Pvalb_9':'Pvalb',
    'SMC':'SMC',
    'Sncg_1':'Sncg','Sncg_2':'Sncg',
    'Sst_1':'Sst','Sst_2':'Sst','Sst_3':'Sst','Sst_4':'Sst','Sst_5':'Sst','Sst_6':'Sst','Sst_7':'Sst','Sst_8':'Sst','Sst_Chodl':'Sst',
    'VLMC':'VLMC',
    'Vip_1':'Vip','Vip_10':'Vip','Vip_2':'Vip','Vip_3':'Vip','Vip_4':'Vip','Vip_5':'Vip','Vip_6':'Vip','Vip_7':'Vip','Vip_8':'Vip','Vip_9':'Vip',
    #'striatum_1','striatum_2',
    #'unannotated',
    #'ventricle_1','ventricle_2'
}
mop_df['ontology'] = mop_df['ontology'].map(mop_to_10x_ann_map)
mop_df = mop_df.drop(columns=['ontology']).dropna()
mop_df.to_csv('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/outputs/gene_cell/MOp_central_10X_ontology.csv',index=False)

#Calculate z from Lyapunov CLT for each gene in each sample
gb_cols = ['sample_id','gene','ontology']

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


#Write out
mop_agg_df.to_csv('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/outputs/gene_ontology/MOp_central_10X_ontology.csv',index=False)

