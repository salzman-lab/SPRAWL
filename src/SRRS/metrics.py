import shapely.geometry

import pandas as pd
import numpy as np
import collections
import itertools
import random
import math
import os

from . import utils
from . import simulate

def _update_med_ranks(cell):
    #mark this cell as ranked to avoid duplicate work
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


def punctate(cell, num_iterations=1000, num_pairs=4):
    """
    Punctate metric

    Steps of this new method:
    1. Remove genes with 1 spot from consideration (should already be done before this)

    2. For each gene, iteratively select X pairs of points and calculate the average distance between them
        Remember this observed average distance for each gene

    3. For P permutations, swap all gene labels
        Repeat step 2 and remember permuted average distances for each gene-count
        Calculate the quantile of the observed average distance against the background of average distances for the corresponding gene counts

    4. Assign a score of 1 if the observed average distance is less than all permutations
        score of -1 if the observed average distance is larger than all permutations, and a score of 0.5 if it is halfway between

    5. Calculate the empirical variance of scores by converting all the null mean dists into scores

    """

    #Taken from https://stackoverflow.com/questions/22229796/choose-at-random-from-combinations
    def random_combination(iterable, r):
        "Random selection from itertools.combinations(iterable, r)"
        pool = tuple(iterable)
        n = len(pool)
        indices = sorted(random.sample(range(n), r))
        return tuple(pool[i] for i in indices)


    def random_mean_pairs_dist(spots):
        """
        Helper function to choose 'num_pairs' pairs of gene spots and calculate the mean distance for each gene
        Input is an array of spots that it will choose from
        Returns a pd.Dataframe
        """
        d = 0
        for _ in range(num_pairs):
            (x1,y1),(x2,y2) = random_combination(spots,2)
            d += math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))

        return d/num_pairs


    all_genes = np.array([g for z in cell.zslices for g in cell.spot_genes[z]])
    all_spots = np.array([xy for z in cell.zslices for xy in cell.spot_coords[z]])

    # Score each gene against the null of the same number of spots
    scores = {
        'gene':[],
        'spot_count':[],
        'score':[],
        'var':[],
    }

    for gene,count in cell.gene_counts.items():

        # Calculate the obs mean dist
        gene_spots = all_spots[all_genes == gene]
        obs_dist = random_mean_pairs_dist(gene_spots)

        # Null distribution by gene-label swapping
        null_dists = []
        for i in range(num_iterations):
            spot_inds = np.random.choice(cell.n,count,replace=False)
            null = random_mean_pairs_dist(all_spots[spot_inds])
            null_dists.append(null)

        null_dists = np.array(null_dists)

        obs = sum(null_dists < obs_dist)
        exp = len(null_dists)/2
        score = (exp-obs)/exp
        null_var = np.var([(exp-d)/exp for d in null_dists])

        scores['gene'].append(gene)
        scores['spot_count'].append(count)
        scores['score'].append(score)
        scores['var'].append(null_var)

    return pd.DataFrame(scores)



def central(cell):
    """
    Central metric

    Not exactly the opposite of the peripheral metric
    Ranks distance from cell centroid
    """
    if cell.ranked:
        return cell

    cell = _central_dist_and_rank(cell)
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
        cell.spot_values[zslice] = min_periph_dists[start_i:end_i]
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
        angle = np.arccos(np.clip(np.dot(norm_vec, mean_gene_vecs[gene]), -1.0, 1.0)) #taken from: https://stackoverflow.com/questions/2827393
        angle_residuals.append(angle)


    #Rank spots by angle residuals
    spot_ranks = np.array(angle_residuals).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

    #save the ranks back to the cell object by z-slice
    start_i = 0
    for zslice in cell.zslices:
        end_i = cell.n_per_z[zslice]+start_i
        cell.spot_ranks[zslice] = spot_ranks[start_i:end_i]
        cell.spot_values[zslice] = angle_residuals[start_i:end_i]
        start_i = end_i

    return cell


def _punctate_dist_and_rank(cell):
    """
    Helper function to calculate punctate ranks
    """
    #gather the spot coords and genes
    spot_genes = []
    spot_coords = []
    for zslice in cell.zslices:
        z_spot_coords = cell.spot_coords[zslice]
        z_spot_genes = cell.spot_genes[zslice]

        spot_coords.extend(z_spot_coords)
        spot_genes.extend(z_spot_genes)

    spot_genes = np.array(spot_genes)
    spot_coords = np.array(spot_coords)

    #calculate gene centroids
    gene_centroids = {}
    for gene in cell.genes:
        gene_inds = spot_genes == gene
        gene_spots = spot_coords[gene_inds]
        gene_centroids[gene] = np.mean(gene_spots,axis=0)

    #calculate per-spot distance to its gene centroid
    spot_dists = []
    for gene,vec in zip(spot_genes,spot_coords):
        cx,cy = gene_centroids[gene]
        x,y = vec
        spot_dists.append((cx-x)*(cx-x)+(cy-y)*(cy-y))

    #Rank spots by angle residuals
    spot_ranks = np.array(spot_dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

    #save the ranks back to the cell object by z-slice
    start_i = 0
    for zslice in cell.zslices:
        end_i = cell.n_per_z[zslice]+start_i
        cell.spot_ranks[zslice] = spot_ranks[start_i:end_i]
        cell.spot_values[zslice] = spot_dists[start_i:end_i]
        start_i = end_i

    return cell


def _central_dist_and_rank(cell):
    """
    Helper function to calculate central ranks
    """

    #calculate distance of spot coords to cell centroid
    spot_genes = []
    spot_dists = []
    for zslice in cell.zslices:
        z_spot_coords = cell.spot_coords[zslice]
        z_spot_genes = cell.spot_genes[zslice]

        #Euclidean distance to slice centroid
        slice_centroid = np.mean(cell.boundaries[zslice], axis=0)
        dists = np.sum((z_spot_coords-slice_centroid)**2, axis=1)

        spot_genes.extend(z_spot_genes)
        spot_dists.extend(dists)

    #Rank the spots
    spot_genes = np.array(spot_genes)
    spot_ranks = np.array(spot_dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

    #save the ranks back to the cell object by z-slice
    start_i = 0
    for zslice in cell.zslices:
        end_i = cell.n_per_z[zslice]+start_i
        cell.spot_ranks[zslice] = spot_ranks[start_i:end_i]
        cell.spot_values[zslice] = spot_dists[start_i:end_i]
        start_i = end_i

    return cell



