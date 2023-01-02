import shapely.geometry

import multiprocessing as mp
import functools

import pandas as pd
import numpy as np
import collections
import itertools
import random
import math
import time
import sys
import os

from . import utils

########################
#   Metrics guidelines #
########################
# Each metric function needs to take in a list of Cell object (from cell.py)
# and output a table of Cell ids with gene_scores and gene_vars

# The metric function calculates the per-gene score for all genes in the cell
#   e.g. the periphery ranking, this will be based on the minimum distance of each spot to the periph
#   e.g. the radial ranking this will be based on the min angle dist of each spot to the gene radial center

def _test(cells, **kwargs):
    """
    Test metric
    Returns a table with all scores of 0 and variances of 1
    """
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

    for cell in cells:
        for gene in cell.genes:
            data['metric'].append('test')
            data['cell_id'].append(cell.cell_id)
            data['annotation'].append(cell.annotation)
            data['num_spots'].append(cell.n)
            data['gene'].append(gene)
            data['num_gene_spots'].append(cell.gene_counts[gene])
            data['score'].append(0)
            data['variance'].append(1)

    return pd.DataFrame(data)


def _peripheral(cell, ret_spot_ranks=False):
    """
    Helper peripheral function gets called by peripheral() for multiprocessing of each cell
    """

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

    #score the cell
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

    #mostly for debugging
    if ret_spot_ranks:
        return spot_genes,spot_ranks


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


def peripheral(cells, **kwargs):
    """
    Peripheral metric

    Steps of this method:
    1. Calculate theoretical gene vars

    2. Measure minimum distance of each spot to the cell-boundary in its z-slice

    3. Rank all spots over all z-slices

    4. Calculate the median gene rank and convert to a score

    5. Return a score dataframe

    """

    #handle kwargs
    processes = kwargs.get('processes', 2)

    #calculate the theoretical gene/cell variances (multiprocessing)
    cells = utils._iter_vars(cells, processes)

    with mp.Pool(processes=processes) as p:
        score_df = pd.concat(p.imap_unordered(_peripheral, cells), ignore_index=True)

    return score_df




def _radial(cell, num_iterations, num_pairs):
    """
    Helper radial function gets called by radial() for multiprocessing of each cell
    """
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

    cell_centroid = np.mean(np.vstack(list(cell.boundaries.values())),axis=0)

    cell = cell.filter_genes_by_count(min_gene_spots=2)

    all_genes = np.array([g for z in cell.zslices for g in cell.spot_genes[z]])
    all_spots = np.array([xy for z in cell.zslices for xy in cell.spot_coords[z]])

    #theoretical variance depends on just the number of iterations
    #Let Y ~ Discrete Uniform from 0 to n, where n is the number of iterations
    #Y represents the number of null permutations that are less than the obs distance and ranges from 0 to 'all'
    #But our score is X = (Y-n/2)/(n/2) because we wanted to center it at 0 and scale it to have values between -1 and 1

    #E[Y] = n/2 from definition of discrete uniform that ranges from 0 to n
    #Var[Y] = n(n+2)/12 

    #E[X] = (2/n)(E[Y]-n/2) = (2/n)(n/2-n/2) = 0
    #Var[X] = (4/n^2)Var[Y] #since Var(aX+b) = (a^2)Var[X]
    #Var[X] = (4/n^2)(n(n+2)/12)
    #Var[X] = (1/n^2)(n(n+2)/3)
    #Var[X] = (1/n)((n+2)/3)
    #Var[X] = (n+2)/3n

    #Also as n --> inf, Var[X] --> 1/3
    #Var[X] = (1+2/n)/3 --> 1/3

    var = (num_iterations+2)/(3*num_iterations)

    pre_calc_nulls = {}

    for gene,count in cell.gene_counts.items():

        # Calculate the obs mean dist
        gene_spots = all_spots[all_genes == gene]
        obs_dist = utils.random_mean_pairs_angle(gene_spots, cell_centroid, num_pairs)

        # Null distribution by gene-label swapping
        if count in pre_calc_nulls:
            null_dists = pre_calc_nulls[count]

        else:
            null_dists = []

            #This is the time-consuming portion, having to make thousands of permutations for the null
            start = time.time()
            for i in range(num_iterations):
                spot_inds = np.random.choice(cell.n,count,replace=False)
                null = utils.random_mean_pairs_angle(all_spots[spot_inds], cell_centroid, num_pairs)
                null_dists.append(null)

            null_dists = np.array(null_dists)
            pre_calc_nulls[count] = null_dists

        obs = sum(null_dists < obs_dist)
        exp = num_iterations/2
        score = (exp-obs)/exp

        data['metric'].append('radial')
        data['cell_id'].append(cell.cell_id)
        data['annotation'].append(cell.annotation)
        data['num_spots'].append(cell.n)
        data['gene'].append(gene)
        data['num_gene_spots'].append(cell.gene_counts[gene])
        data['score'].append(score)
        data['variance'].append(var)


    return pd.DataFrame(data)



def radial(cells, **kwargs):
    """
    Radial metric

    Steps of this method:
    1. Remove genes with 1 spot from consideration

    2. For each gene, iteratively select X pairs of points and calculate the average angle between them
       the angle being measured formed by (x1,y1) --> (cx,cy) --> (x2,y2) where (cx,cy) is the cell centroid

    3. For `num_iterations` permutations, swap all gene labels
       Repeat step 2 and remember permuted average angle for each gene-count
       Calculate the quantile of the observed average distance against the background of average distances for the corresponding gene counts

    4. Assign a score of 1 if the observed average distance is less than all permutations
       score of -1 if the observed average distance is larger than all permutations, and a score of 0.5 if it is halfway between

    5. Calculate the empirical variance of scores by converting all the null mean angles into scores

    """
    #handle kwargs
    processes = kwargs.get('processes', 2)
    num_iterations = kwargs.get('num_iterations', 1000)
    num_pairs = kwargs.get('num_pairs', 4)

    f = functools.partial(_radial, num_iterations=num_iterations, num_pairs=num_pairs)
    
    with mp.Pool(processes=processes) as p:
        score_df = pd.DataFrame()
        for i,cell_df in enumerate(p.imap_unordered(f, cells)):
            score_df = pd.concat((score_df, cell_df), ignore_index=True)

    return score_df


def _punctate(cell, num_iterations, num_pairs):
    """
    Helper punctate function gets called by punctate() for multiprocessing of each cell
    """
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

    cell = cell.filter_genes_by_count(min_gene_spots=2)

    all_genes = np.array([g for z in cell.zslices for g in cell.spot_genes[z]])
    all_spots = np.array([xy for z in cell.zslices for xy in cell.spot_coords[z]])

    #theoretical variance depends on just the number of iterations
    #Let Y ~ Discrete Uniform from 0 to n, where n is the number of iterations
    #Y represents the number of null permutations that are less than the obs distance and ranges from 0 to 'all'
    #But our score is X = (Y-n/2)/(n/2) because we wanted to center it at 0 and scale it to have values between -1 and 1

    #E[Y] = n/2 from definition of discrete uniform that ranges from 0 to n
    #Var[Y] = n(n+2)/12 

    #E[X] = (2/n)(E[Y]-n/2) = (2/n)(n/2-n/2) = 0
    #Var[X] = (4/n^2)Var[Y] #since Var(aX+b) = (a^2)Var[X]
    #Var[X] = (4/n^2)(n(n+2)/12)
    #Var[X] = (1/n^2)(n(n+2)/3)
    #Var[X] = (1/n)((n+2)/3)
    #Var[X] = (n+2)/3n
    var = (num_iterations+2)/(3*num_iterations)

    pre_calc_nulls = {}
    for gene,count in cell.gene_counts.items():

        # Calculate the obs mean dist
        gene_spots = all_spots[all_genes == gene]
        obs_dist = utils.random_mean_pairs_dist(gene_spots, num_pairs)

        # Null distribution by gene-label swapping
        if count in pre_calc_nulls:
            null_dists = pre_calc_nulls[count]
        else:
            null_dists = []
            for i in range(num_iterations):
                spot_inds = np.random.choice(cell.n,count,replace=False)
                null = utils.random_mean_pairs_dist(all_spots[spot_inds], num_pairs)
                null_dists.append(null)

            null_dists = np.array(null_dists)
            pre_calc_nulls[count] = null_dists

        obs = sum(null_dists < obs_dist)
        exp = len(null_dists)/2
        score = (exp-obs)/exp


        data['metric'].append('puncta')
        data['cell_id'].append(cell.cell_id)
        data['annotation'].append(cell.annotation)
        data['num_spots'].append(cell.n)
        data['gene'].append(gene)
        data['num_gene_spots'].append(cell.gene_counts[gene])
        data['score'].append(score)
        data['variance'].append(var)

    return pd.DataFrame(data)


def punctate(cells, **kwargs):
    """
    Punctate metric

    Uses multiprocessing to separately work on each cell
    Steps of this new method:
    1. Remove genes with 1 spot from consideration

    2. For each gene, iteratively select X pairs of points and calculate the average distance between them
        Remember this observed average distance for each gene

    3. For `num_iterations` permutations, swap all gene labels
        Repeat step 2 and remember permuted average distances for each gene-count
        Calculate the quantile of the observed average distance against the background of average distances for the corresponding gene counts

    4. Assign a score of 1 if the observed average distance is less than all permutations
        score of -1 if the observed average distance is larger than all permutations, and a score of 0.5 if it is halfway between

    5. Calculate the empirical variance of scores by converting all the null mean dists into scores

    """
    #handle kwargs
    processes = kwargs.get('processes', 2)
    num_iterations = kwargs.get('num_iterations', 1000)
    num_pairs = kwargs.get('num_pairs', 4)

    f = functools.partial(_punctate, num_iterations=num_iterations, num_pairs=num_pairs)
    
    with mp.Pool(processes=processes) as p:
        score_df = pd.concat(p.imap_unordered(f, cells), ignore_index=True)

    return score_df


def _central(cell):
    """
    Helper central function gets called by central() for multiprocessing of each cell
    """
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

    #score the genes
    exp_med = (cell.n+1)/2
    for gene in cell.genes:
        gene_ranks = spot_ranks[spot_genes == gene]
        obs_med = np.median(gene_ranks)
        score = (exp_med-obs_med)/(exp_med-1)

        data['metric'].append('central')
        data['cell_id'].append(cell.cell_id)
        data['annotation'].append(cell.annotation)
        data['num_spots'].append(cell.n)
        data['gene'].append(gene)
        data['num_gene_spots'].append(cell.gene_counts[gene])
        data['score'].append(score)
        data['variance'].append(cell.gene_vars[gene])

    return pd.DataFrame(data)



def central(cells, **kwargs):
    """
    Central metric

    Not exactly the opposite of the peripheral metric
    Ranks distance from cell centroid
    """
    processes = kwargs.get('processes', 2)

    #calculate the theoretical gene/cell variances (multiprocessing)
    cells = utils._iter_vars(cells, processes)

    with mp.Pool(processes=processes) as p:
        score_df = pd.concat(p.imap_unordered(_central, cells), ignore_index=True)

    return score_df



