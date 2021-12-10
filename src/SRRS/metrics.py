import shapely.geometry
from scipy import stats

import pandas as pd
import numpy as np
import collections
import h5py
import os

import operator as op
from functools import reduce
import math
import time

from . import utils


########################
#   Metrics guidelines #
########################
# Each metric needs to take input in the form of an hdf5_path and a cell_id
# The metric function calculates the per-gene score for all genes in the cell
#   e.g. the periphery ranking, this will be based on the minimum distance of each spot to the periph
#   e.g. the radial ranking this will be based on the min angle dist of each spot to the gene radial center

#Return value is expected to be a pandas dataframe where each row is a unique gene with the following required columns
#     'metric' peripheral/polar etc
#     'cell_id' the cell id. will be the same for every row
#     'annotation' The cell-type annotation information. same for all rows
#     'num_spots' The number of spots in the whole cell. same for all rows
#     'gene' The gene of the row
#     'num_gene_spots' The number of spots of this spots gene in the cell
#     'median_rank' the observed median rank
#     'score' value ranging from -1 to 1
#     'variance' theoretical variance under the null

#extra columns shouldn't be a problem

def _test(hdf5_path, cell_id):
    """
    Test metric
    Returns a test df with all columns as expected except
        metric - test
        score - 0
    """

    with h5py.File(hdf5_path,'r') as f:
        cell = f['cells'][cell_id]
        gene_counts = collections.defaultdict(int)

        for z_ind in cell.attrs.get('zslices'):
            z_counts = collections.Counter(cell['spot_genes'][z_ind])
            gene_counts.update(z_counts)

        #drop a spot for genes with an even number of spots
        n = sum(gene_counts.values())
        ms = [m if m%2 == 1 else m-1 for m in gene_counts.values()]
        vs = [utils.calc_var(m,n) for m in ms]

        df = pd.DataFrame({
            'metric':'test',
            'cell_id':cell_id,
            'annotation':cell.attrs['annotation'],
            'num_spots':n,
            'gene':gene_counts.keys(),
            'num_gene_spots':ms,
            'median_rank':(n+1)/2,
            'score':0,
            'variance':vs,
        })

    return df


def peripheral(hdf5_path, cell_id):
    """
    Peripheral metric
    """
    min_periph_dists = []
    min_spot_genes = []

    with h5py.File(hdf5_path,'r') as f:
        cell = f['cells'][cell_id]
        annotation = cell.attrs.get('annotation')

        for zslice in cell.attrs.get('zslices'):

            #Calculate dists of each spot to periphery
            boundary = cell['boundaries'][zslice]
            spot_coords = cell['spot_coords'][zslice]
            spot_genes = cell['spot_genes'][zslice]

            poly = shapely.geometry.Polygon(boundary)
            for p,gene in zip(spot_coords,spot_genes):
                dist = poly.boundary.distance(shapely.geometry.Point(p))
                min_periph_dists.append(dist)
                min_spot_genes.append(gene)

    #Rank the spots
    min_spot_genes = np.array(min_spot_genes)
    spot_ranks = np.array(min_periph_dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0
    tot_spots = len(spot_ranks)

    #Iterate through unique genes to get per-gene score
    genes = np.unique(min_spot_genes)
    num_genes_per_spot = []
    gene_scores = []
    mn_vars = []
    obs_medians = []
    for gene in genes:
        gene_inds = min_spot_genes == gene
        num_spots = int(sum(gene_inds))

        gene_ranks = spot_ranks[gene_inds]

        #Drop one spot randomly (uniformly) if an even number
        #approximation. reason is calculating variance for even number of order of magnitude slower
        if num_spots%2 == 0:
            num_spots -= 1
            gene_ranks = np.random.choice(gene_ranks,num_spots,replace=False)

        num_genes_per_spot.append(num_spots)
        obs_med = np.median(gene_ranks)
        score = utils.score(obs_med, tot_spots)
        var = utils.calc_var(num_spots,tot_spots)

        gene_scores.append(score)
        mn_vars.append(var)
        obs_medians.append(obs_med)

    df = pd.DataFrame({
        'metric':'peripheral',
        'cell_id':cell_id,
        'annotation':annotation,
        'num_spots':tot_spots,
        'gene':[g.decode() for g in genes],
        'num_gene_spots':num_genes_per_spot,
        'median_rank':obs_medians,
        'score':gene_scores,
        'variance':mn_vars,
    })

    return df




