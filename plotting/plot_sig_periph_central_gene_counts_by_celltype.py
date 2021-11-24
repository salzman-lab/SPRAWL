#!/bin/bash
#############################
# File Name : plot_sig_periph_central_gene_counts_by_celltype.py
#
# Purpose : [???]
#
# Creation Date : 24-11-2021
#
# Last Modified : Wed 24 Nov 2021 11:15:05 AM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

#Read in all the info at the gene/cell-type level
f_paths = glob.glob('/oak/stanford/groups/horence/rob/isoform_localizations/sub_cell_patterns/outputs/20211104_merfish_periph_zscores_one_drop_evens_converted_celltypes/m*final_output.csv')
merf_df = pd.DataFrame()

for f_path in f_paths:
    sample_name = os.path.basename(f_path).split('_')[1]
    sub_df = pd.read_csv(f_path)
    sub_df['sample_name'] = sample_name
    merf_df = pd.concat((merf_df, sub_df), ignore_index=True)

merf_df['direction'] = np.where(merf_df['z_score'].le(0),'central','peripheral')
merf_df.head()


sig_rows = merf_df['bh_corrected_two_sided_p'].le(0.05)
pos_rows = merf_df['z_score'].ge(0)

merf_df.loc[~sig_rows, 'Category'] = 'Insignificant'
merf_df.loc[sig_rows & pos_rows, 'Category'] = 'Peripheral'
merf_df.loc[sig_rows & (~pos_rows), 'Category'] = 'Central'

scatter_df = merf_df.groupby(['sample_name','annotation','Category']).size().unstack().reset_index()


g = sns.relplot(
    data=scatter_df, 
    x='Peripheral', 
    y='Central',
    hue='annotation', 
    palette='tab20',
    col='sample_name',
    height=2,
    col_wrap=3,
)

for ax in g.axes_dict.values():
    ax.axline((0, 0), slope=1, c=".2", ls="--", zorder=0, color='grey')
    ax.set_xlim(0,50)
    ax.set_ylim(0,50)
    
    
g.fig.subplots_adjust(top=0.9)
plt.suptitle('Number of significant peripheral/central genes per cell-type per sample')
#plt.show() #change to a plt.savefig(...)
plt.savefig('scatter_sig_gene_counts_plot.png')
plt.close()

