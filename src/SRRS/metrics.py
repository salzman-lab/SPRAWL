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
# Each metric function needs to take in a Cell object (from cell.py)
# and output the same Cell with gene_med_ranks dictionary filled in for each gene
# The metric function calculates the per-gene score for all genes in the cell
#   e.g. the periphery ranking, this will be based on the minimum distance of each spot to the periph
#   e.g. the radial ranking this will be based on the min angle dist of each spot to the gene radial center

def _test(cell):
    """
    Test metric
    Returns a cell with all med ranks equal to (cell.n+1)/2
        score - 0
    """
    for g in cell.genes:
        cell.gene_med_ranks[g] = (cell.n+1)/2

    return cell


def peripheral(cell):
    """
    Peripheral metric
    """
    #Use helper function to calculate distances and ranks
    min_spot_genes,spot_ranks = _peripheral_dist_and_rank(cell)

    #Iterate through unique genes to assign per-gene scores
    for gene in cell.genes:
        gene_inds = min_spot_genes == gene
        gene_ranks = spot_ranks[gene_inds]
        cell.gene_med_ranks[gene] = np.median(gene_ranks)

    return cell


def _peripheral_dist_and_rank(cell):
    """
    Helper function to calculate peripheral ranks
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

    return min_spot_genes,spot_ranks


