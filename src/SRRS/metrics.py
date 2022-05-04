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

########################
#   Metrics guidelines #
########################
# Each metric function needs to take in a list of Cell object (from cell.py)
# and output the a table of Cell ids with gene_scores and gene_vars

# The metric function calculates the per-gene score for all genes in the cell
#   e.g. the periphery ranking, this will be based on the minimum distance of each spot to the periph
#   e.g. the radial ranking this will be based on the min angle dist of each spot to the gene radial center

def _test(cells):
    """
    Test metric
    Returns a cell with all med ranks equal to (cell.n+1)/2
        score - 0
    """
    for g in cell.genes:
        cell.gene_scores[g] = 0
        cell.gene_vars[g] = 1

    return cell


def peripheral(cells, processes=2):
    """
    Peripheral metric
    """

    #calculate the theoretical gene/cell variances (multiprocessing)
    cells = utils._iter_vars(cells, processes)

    data = {
        'metric':[],
        'cell_id':[],
        'annotation':[],
        'num_spots':[],
        'gene':[],
        'num_gene_spots':[],
        'score':[],
        'variance':[],
    }

    #score the cells (NOTE! make this multiprocessing)
    for cell in cells:
        periph_dists = []
        spot_genes = []

        for zslice in cell.zslices:

            #Calculate dists of each spot to periphery
            z_boundary = cell.boundaries[zslice]
            z_spot_coords = cell.spot_coords[zslice]
            z_spot_genes = cell.spot_genes[zslice]

            poly = shapely.geometry.Polygon(z_boundary)
            for xy,gene in zip(z_spot_coords,z_spot_genes):
                dist = poly.boundary.distance(shapely.geometry.Point(xy))
                periph_dists.append(dist)
                spot_genes.append(gene)

        #Rank the spots
        spot_genes = np.array(spot_genes)
        spot_ranks = np.array(periph_dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

        #score the genes
        exp_med = (cell.n+1)/2
        for gene in cell.genes:
            gene_ranks = spot_ranks[spot_genes == gene]
            obs_med = np.median(gene_ranks)
            score = (exp_med-obs_med)/(exp_med-1)

            data['metric'].append('periph')
            data['cell_id'].append(cell.cell_id)
            data['annotation'].append(cell.annotation)
            data['num_spots'].append(cell.n)
            data['gene'].append(gene)
            data['num_gene_spots'].append(cell.gene_counts[gene])
            data['score'].append(score)
            data['variance'].append(cell.gene_vars[gene])

    return pd.DataFrame(data)


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

    pre_calc_nulls = {}

    for gene,count in cell.gene_counts.items():

        # Calculate the obs mean dist
        gene_spots = all_spots[all_genes == gene]
        obs_dist = random_mean_pairs_dist(gene_spots)

        # Null distribution by gene-label swapping
        if count in pre_calc_nulls:
            null_dists = pre_calc_nulls[count]

        else:
            null_dists = []
            for i in range(num_iterations):
                spot_inds = np.random.choice(cell.n,count,replace=False)
                null = random_mean_pairs_dist(all_spots[spot_inds])
                null_dists.append(null)

            null_dists = np.array(null_dists)
            pre_calc_nulls[count] = null_dists

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



