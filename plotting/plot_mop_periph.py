#!/bin/bash
#############################
# File Name : plot_mop_periph.py
#
# Purpose : [???]
#
# Creation Date : 23-11-2021
#
# Last Modified : Tue 23 Nov 2021 01:07:11 PM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import glob
import sys
import os

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

#Params
gene = sys.argv[1]
min_samps_of_ann = 1


#Read in all the info at the gene/cell-type level
f_paths = glob.glob('/oak/stanford/groups/horence/rob/isoform_localizations/sub_cell_patterns/outputs/20211104_merfish_periph_zscores_one_drop_evens_converted_celltypes/m*final_output.csv')
merf_df = pd.DataFrame()

for f_path in f_paths:
    sample_name = os.path.basename(f_path).split('_')[1]
    sub_df = pd.read_csv(f_path)
    sub_df['sample_name'] = sample_name
    merf_df = pd.concat((merf_df, sub_df), ignore_index=True)

merf_df['direction'] = np.where(merf_df['z_score'].le(0),'central','peripheral')

#Read in all the info at the cell level
f_paths = glob.glob('/oak/stanford/groups/horence/rob/isoform_localizations/sub_cell_patterns/outputs/20211104_merfish_periph_zscores_one_drop_evens_converted_celltypes/m*_filtered_cells.csv')
cells_df = pd.DataFrame()

for f_path in f_paths:
    sample_name = os.path.basename(f_path).split('_')[1]
    sub_df = pd.read_csv(f_path)
    sub_df['sample_name'] = sample_name
    cells_df = pd.concat((cells_df, sub_df), ignore_index=True)




gene_ann_samps = (
    merf_df[merf_df['gene'].eq(gene)][
        ['gene','annotation','sample_name','bh_corrected_two_sided_p','num_cells']
    ]
)

gene_ann_samps = gene_ann_samps.groupby('annotation').filter(lambda g: len(g) >= min_samps_of_ann)
gene_ann_samps['Mouse'] = np.where(gene_ann_samps['sample_name'].str.startswith('M1'),'Mouse1','Mouse2')
ann_order = sorted(gene_ann_samps['annotation'].unique())

#Iterate through the two mice separately
fig,axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(10,5))

for i,(mouse,mouse_gene_ann_samps) in enumerate(gene_ann_samps.groupby('Mouse')):
    print(mouse)
    
    plot_df = mouse_gene_ann_samps.merge(
        cells_df,
        on=['gene','annotation','sample_name'],
    ).sort_values(['annotation','sample_name'])


    sig_box = plot_df.drop_duplicates(['annotation','sample_name'])['bh_corrected_two_sided_p'].le(0.05)
    sig_box = list(sig_box)


    frac_bins = np.arange(0,1,0.1)
    gene_fracs = plot_df['num_gene_spots'].div(plot_df['num_spots'])
    plot_df['frac_bin'] = pd.cut(gene_fracs, frac_bins)
    
    ##################
    #                #
    #     Scatters   #
    #                #
    ##################
    gfs = sorted(plot_df['frac_bin'].unique())
    gf_colors = sns.color_palette("rocket",len(gfs))

    legend_elems = []

    for gf,gf_color in zip(gfs,gf_colors):

        sub_plot_df = plot_df[plot_df['frac_bin'].eq(gf)]

        strip_ax = sns.stripplot(
            x = 'annotation',
            y = 'gene_score',
            hue = 'sample_name',
            palette = [gf_color],
            alpha = 0.5,
            dodge = True,
            order = ann_order,
            data = sub_plot_df,
            ax = axs[i],
        )

        legend_elems.append(Patch(facecolor=gf_color,label=gf))

    ##################
    #                #
    #     Boxplots   #
    #                #
    ##################
    box_ax = sns.boxplot(
        x = 'annotation',
        y = 'gene_score',
        hue = 'sample_name',
        color = 'grey',
        fliersize = 0,
        order = ann_order,
        data = plot_df,
        ax = axs[i],
    )

    #Adding stars above boxplots if significant
    lines = box_ax.get_lines()
    boxes = [c for c in box_ax.get_children() if type(c).__name__ == 'PathPatch']
    lines_per_box = int(len(lines) / len(boxes))
    for box_num,median in enumerate(lines[4:len(lines):lines_per_box]):
        t = '*' if sig_box[box_num] else ''
        x, y = (data.mean() for data in median.get_data())
        text = box_ax.text(x, 1.1, t, ha='center', va='center', color='black', fontweight='bold')

    #Subplot controlling look
    axs[i].set_ylabel(mouse)
    axs[i].legend().set_visible(False)
    axs[i].set_xlabel('')
    axs[i].set_ylim(-1.2,1.2) #little bit larger so asterixis show
    axs[i].set_yticks([-1.0,-0.5,0.0,0.5,1.0])
    axs[i].axhline(0,color='grey',linestyle='dotted')

##################
#                #
#      Legend    #
#                #
##################
# create legend just for boxplot to avoid duplicates
l = plt.legend(
    title = 'Gene fraction per cell',
    handles=legend_elems,
    bbox_to_anchor=(1.05, 1),
    loc=2,
    borderaxespad=0.,
)

fig.text(-0.03, 0.5, 'Effect size (positive is peripheral)', va='center', rotation='vertical')
plt.suptitle(gene)
plt.tight_layout()

#plt.show()
plt.savefig('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/plotting/no_min_sample_boxplots/{}.png'.format(gene))
plt.close()

