import shapely.geometry
from scipy import stats

import pandas as pd
import numpy as np
import collections
import h5py
import os

import operator as op
from functools import reduce
import logging
import math
import time

from . import utils


########################
#   Metrics guidelines #
########################
# Each metric function needs to have the following method signature
#   def metric_name(cell, var_mem={})

# First argument is a Cell object (from cell.py)
# Second optional argument var_mem is a dict to cache theoretical vars
#
#
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

def _test(cell, var_mem={}):
    """
    Test metric
    Returns a test df with all columns as expected except
        metric - test
        score - 0
    """

    n = cell.n
    ms = cell.gene_counts.values()

    vs = []
    for m in ms:
        if (m,n) not in var_mem:
            var_mem[(m,n)] = utils.calc_var(m,n)
        vs.append(var_mem[(m,n)])

    df = pd.DataFrame({
        'metric':'test',
        'cell_id':cell.cell_id,
        'annotation':cell.annotation,
        'num_spots':cell.n,
        'gene':[g.decode() for g in cell.gene_counts.keys()],
        'num_gene_spots':ms,
        'median_rank':(n+1)/2,
        'score':0,
        'variance':vs,
    })

    return df


def peripheral(cell, var_mem={}):
    """
    Peripheral metric
    """
    min_periph_dists = []
    min_spot_genes = []

    for zslice in cell.zslices:

        #Calculate dists of each spot to periphery
        boundary = cell.boundaries[zslice]
        spot_coords = cell.spot_coords[zslice]
        spot_genes = cell.spot_genes[zslice]

        poly = shapely.geometry.Polygon(boundary)
        for p,gene in zip(spot_coords,spot_genes):
            dist = poly.boundary.distance(shapely.geometry.Point(p))
            min_periph_dists.append(dist)
            min_spot_genes.append(gene)

    #Rank the spots
    min_spot_genes = np.array(min_spot_genes)
    spot_ranks = np.array(min_periph_dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0
    n = len(spot_ranks)

    #Iterate through unique genes to get per-gene score
    genes = np.unique(min_spot_genes)
    num_genes_per_spot = []
    gene_scores = []
    mn_vars = []
    obs_medians = []
    for gene in genes:
        gene_inds = min_spot_genes == gene
        m = int(sum(gene_inds))

        gene_ranks = spot_ranks[gene_inds]

        num_genes_per_spot.append(m)
        obs_med = np.median(gene_ranks)
        score = utils.score(obs_med, n)

        #calculate var if not cached
        if (m,n) not in var_mem:
            var_mem[(m,n)] = utils.calc_var(m,n)
        var = var_mem[(m,n)]

        gene_scores.append(score)
        mn_vars.append(var)
        obs_medians.append(obs_med)

    df = pd.DataFrame({
        'metric':'peripheral',
        'cell_id':cell.cell_id,
        'annotation':cell.annotation,
        'num_spots':n,
        'gene':[g.decode() for g in genes],
        'num_gene_spots':num_genes_per_spot,
        'median_rank':obs_medians,
        'score':gene_scores,
        'variance':mn_vars,
    })

    return df




