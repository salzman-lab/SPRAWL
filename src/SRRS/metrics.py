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

def _update_cell(cell,genes,ranks):
    #save the ranks back to the cell object by z-slice
    start_i = 0
    for zslice in cell.zslices:
        end_i = cell.n_per_z[zslice]+start_i
        cell.spot_ranks[zslice] = ranks[start_i:end_i]
        start_i = end_i

    cell.ranked = True

    #Iterate through unique genes to assign per-gene scores
    for gene in cell.genes:
        gene_inds = genes == gene
        gene_ranks = ranks[gene_inds]
        cell.gene_med_ranks[gene] = np.median(gene_ranks)

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

    genes,ranks = _peripheral_dist_and_rank(cell)
    cell = _update_cell(cell,genes,ranks)
    return cell


def radial(cell):
    """
    Radial metric
    """
    if cell.ranked:
        return cell

    genes,ranks = _radial_dist_and_rank(cell)
    cell = _update_cell(cell,genes,ranks)
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



def _radial_dist_and_rank(cell):
    """
    Helper function to calculate radial ranks
    """

    #calculate the spot angles
    spot_angles = []
    spot_genes = []

    for zslice in cell.zslices:
        z_boundary = cell.boundaries[zslice]
        z_spot_coords = cell.spot_coords[zslice]
        z_spot_genes = cell.spot_genes[zslice]

        z_centroid = np.mean(z_boundary,axis=0)
        centered_spots = z_spot_coords-z_centroid

        x = centered_spots[:,0]
        y = centered_spots[:,1]
        z_horiz_angs = np.abs(np.arctan2(y,x))

        spot_angles.extend(z_horiz_angs)
        spot_genes.append(z_spot_genes)


    #Calculate median gene angles and residuals from median
    angle_residuals = spot_angles.copy()
    for gene in cell.genes:
        gene_inds = spot_genes == gene
        gene_angles = spot_angles[gene_inds]
        med_gene_angle = np.median(gene_angles)
        angle_residuals[gene_inds] = np.abs(gene_angles-med_gene_angle)

    #Rank spots by angle residuals
    spot_ranks = np.array(angle_residuals).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

    #save the ranks back to the cell object by z-slice
    for zslice in cell.zslices:
        pass

    return spot_genes,spot_ranks


