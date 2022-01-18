import shapely.geometry
from scipy import stats

import pandas as pd
import numpy as np
import collections
import h5py
import os

import operator as op
from functools import reduce
import collections
import logging
import math
import time

from . import utils

def _update_med_ranks(cell):
    #This cell must have been ranked
    cell.ranked = True

    #Pull out the gene ranks
    gene_ranks = collections.defaultdict(list)
    for z in cell.zslices:
        for gene,rank in zip(cell.spot_genes[z],cell.spot_ranks[z]):
            gene_ranks[gene].append(rank)

    #Iterate through unique genes to assign per-gene scores
    for gene,ranks in gene_ranks.items():
        cell.gene_med_ranks[gene] = np.median(ranks)

    return cell


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
    if cell.ranked:
        return cell

    cell = _peripheral_dist_and_rank(cell)
    cell = _update_med_ranks(cell)
    return cell


def radial(cell):
    """
    Radial metric
    """
    if cell.ranked:
        return cell

    cell = _radial_dist_and_rank(cell)
    cell = _update_med_ranks(cell)
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

    #save the ranks back to the cell object by z-slice
    start_i = 0
    for zslice in cell.zslices:
        end_i = cell.n_per_z[zslice]+start_i
        cell.spot_ranks[zslice] = spot_ranks[start_i:end_i]
        start_i = end_i

    return cell



def _radial_dist_and_rank(cell):
    """
    Helper function to calculate radial ranks
    """

    #calculate the cell (x,y) centroid over all z-slices
    cell_centroid = np.mean(np.vstack(list(cell.boundaries.values())),axis=0)

    #normalize spot coords to cell centroid
    spot_genes = []
    spot_coords = []
    for zslice in cell.zslices:
        z_spot_coords = cell.spot_coords[zslice]
        z_spot_genes = cell.spot_genes[zslice]

        centered_spots = z_spot_coords-cell_centroid

        spot_coords.extend(centered_spots)
        spot_genes.extend(z_spot_genes)

    spot_genes = np.array(spot_genes)
    spot_coords = np.array(spot_coords)

    #calculate mean gene vecs
    mean_gene_vecs = {}
    for gene in cell.genes:
        gene_inds = spot_genes == gene
        gene_spots = spot_coords[gene_inds]
        mean_vec = np.mean(gene_spots,axis=0)
        mean_gene_vecs[gene] = mean_vec/np.linalg.norm(mean_vec)

    #calculate angle residuals
    angle_residuals = []
    for gene,vec in zip(spot_genes,spot_coords):
        norm_vec = vec/np.linalg.norm(vec)
        dp = np.dot(mean_gene_vecs[gene], vec)
        angle = np.arccos(dp)
        angle_residuals.append(angle)


    #Rank spots by angle residuals
    spot_ranks = np.array(angle_residuals).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

    #save the ranks back to the cell object by z-slice
    start_i = 0
    for zslice in cell.zslices:
        end_i = cell.n_per_z[zslice]+start_i
        cell.spot_ranks[zslice] = spot_ranks[start_i:end_i]
        start_i = end_i

    return cell


