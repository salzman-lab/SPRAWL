#!/bin/bash
#############################
# File Name : create_hdf5.py
#
# Purpose : [???]
#
# Creation Date : 12-11-2021
#
# Last Modified : Mon 22 Nov 2021 11:34:38 AM PST
#
# Created By : Rob Bierman
#
##############################

import h5py
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os

import shapely.geometry

seqfish_stem = '/oak/stanford/groups/horence/rob/SeqFishPlus_data'

############################################
#                                          #
#                                          #
#             Prepare metadata             #
#                                          #
#                                          #
############################################
#Use the header row of the counts file to get the gene names
#Looks like the gene names are 1-indexed since max value is 10000, min is 1 in the RNA_location CSVs
with open(os.path.join(seqfish_stem,'sourcedata','cortex_svz_counts.csv'),'r') as gene_f:
    genes = gene_f.readline()
    gene_ind_to_name = {i+1:g for i,g in enumerate(genes.split(','))}

#Prepare dict to convert from cell_ind to annotation ind
ann_df = pd.read_csv(
    os.path.join(seqfish_stem,'cortex_svz_cell_type_annotations.csv'),
)
cell_id_to_louvain = dict(ann_df[['cell_id','louvain']].values)
cell_id_to_ann = dict(ann_df[['cell_id','annotation']].values)

#Prepare the hdf5 file for write out
f = h5py.File(
    '/oak/stanford/groups/horence/rob/SeqFishPlus_data/seqfish_plus.hdf5',
    mode='w',
)

hdf5_cells = f.create_group('cells')
overall_unique_genes = set()
cell_ids = []


############################################
#                                          #
#                                          #
#         Iterate through cells            #
#                                          #
#                                          #
############################################
rna_loc_paths = glob.glob(os.path.join(seqfish_stem,'sourcedata','RNA_locations','*.csv'))

for i,rna_loc_path in enumerate(rna_loc_paths):
    
    #cell_inds are 1-indexed in the file names, but 0-indexed in the annotation
    #subtracting 1 here so that the indexing is the same
    cell_ind = int(os.path.basename(rna_loc_path).split('_')[-1].split('.')[0])-1
    cell_id = 'cell_{}'.format(cell_ind)
    
    #The table is super wide instead of long, have to transpose it
    df = pd.read_csv(rna_loc_path,header=None).T
    df.columns = ['gene_ind','x','y']
    df['gene_ind'] = df['gene_ind'].astype(int)
    df['gene_name'] = df['gene_ind'].map(gene_ind_to_name)
    
    encoded_genes = list(df['gene_name'].apply(lambda s: s.encode()))
    
    #Calculate the convex hull boundary coordinates
    #Help from https://gis.stackexchange.com/questions/360916/computing-convex-hull-of-points-using-shapely
    rna_spots = shapely.geometry.MultiPoint([
        shapely.geometry.Point(x,y) for x,y in df[['x','y']].values
    ])
    
    boundary = np.array(rna_spots.convex_hull.boundary.xy).T
    
    #Update global genes and cell_id info
    cell_ids.append(cell_id.encode())
    overall_unique_genes = overall_unique_genes.union(encoded_genes)
    
    #Add this info to the hdf5 file
    hdf5_cell = hdf5_cells.create_group(cell_id)
    
    hdf5_cell.attrs['louvain'] = cell_id_to_louvain[cell_ind]
    hdf5_cell.attrs['annotation'] = cell_id_to_ann[cell_ind]
    hdf5_cell.attrs['num_genes'] = df['gene_name'].nunique()
    hdf5_cell.attrs['num_spots'] = len(df)
    hdf5_cell.attrs['zslices'] = ['0'] #just a single z-slice in this dataset

    bounds = hdf5_cell.create_group('boundaries')
    spot_coords = hdf5_cell.create_group('spot_coords')
    spot_genes = hdf5_cell.create_group('spot_genes')
    
    bounds.create_dataset('0', data=boundary)
    spot_coords.create_dataset('0', data=df[['x','y']].values)
    spot_genes.create_dataset('0', data=encoded_genes)

genes = sorted(list(overall_unique_genes))
f.create_dataset('cell_ids',data=cell_ids)
f.create_dataset('genes',data=genes)
f.close()

