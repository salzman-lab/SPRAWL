import scanpy as sc
import pandas as pd
import numpy as np
import os

import seaborn as sns
import matplotlib.pyplot as plt
import time
import h5py
import sys

cbg_paths = [
    ('Liver1','Slice1','/scratch/groups/horence/rob/data/vz_liver_showcase/Liver1Slice1/cell_by_gene.csv'),
    ('Liver1','Slice2','/scratch/groups/horence/rob/data/vz_liver_showcase/Liver1Slice2/cell_by_gene.csv'),
    ('Liver2','Slice1','/scratch/groups/horence/rob/data/vz_liver_showcase/Liver2Slice1/cell_by_gene.csv'),
    ('Liver2','Slice2','/scratch/groups/horence/rob/data/vz_liver_showcase/Liver2Slice2/cell_by_gene.csv'),
]

neighs,pcs = 60,10 #leads to 100 clusters
#neighs,pcs = 40,10 #leads to 100 clusters
#neighs,pcs = 30,20 #leads to 167 clusters
#neighs,pcs = 30,10 #leads to 189 clusters
#neighs,pcs = 15,10 #leads to 207 clusters
#neighs,pcs = 10,10 #leads to 682 clusters
#neighs,pcs = 5,10 #leads to 1928 clusters

df = pd.DataFrame()
samps = []
reps = []

#Read in all the data from all replicates then subset
for sample,rep_num,cbg_path in cbg_paths:
    print(cbg_path)
    cbg_df = pd.read_csv(cbg_path)

    #Subset to the cells with at least 100 spots
    cbg_df = pd.read_csv(cbg_path, index_col=0)
    spots_per_cell = cbg_df.sum(axis=1)
    keep_inds = spots_per_cell >= 100
    cbg_df = cbg_df.loc[keep_inds]

    samps += [sample]*len(cbg_df)
    reps += [rep_num]*len(cbg_df)

    df = pd.concat((df, cbg_df))

print('Size of cell-by-gene table',df.shape)
    
adata = sc.AnnData(df)
adata.obs['sample'] = samps
adata.obs['replicate'] = reps

start = time.time()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print('Normalized',time.time()-start)
sys.stdout.flush()

start = time.time()
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
print('Highly variable genes selected',time.time()-start)
sys.stdout.flush()

start = time.time()
sc.tl.pca(adata, svd_solver='arpack')
print('PCA calculated',time.time()-start)
sys.stdout.flush()

start = time.time()
sc.pp.neighbors(adata, n_neighbors=neighs, n_pcs=pcs)
print('Neighbors calculated',time.time()-start)
sys.stdout.flush()

start = time.time()
sc.tl.leiden(adata)
print('Leiden calculated',time.time()-start)
sys.stdout.flush()

print('Number of clusters',adata.obs['leiden'].nunique())
sys.stdout.flush()

#start = time.time()
#sc.tl.umap(adata)
#sc.pl.umap(adata, color=['leiden','replicate'])
#print('umap calculated and plotted',time.time()-start)

adata.obs.to_csv('vz_liver_cluster_assignments_{}_neighs_{}_pcs.csv'.format(neighs,pcs))


