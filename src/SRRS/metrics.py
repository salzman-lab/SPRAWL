import shapely.geometry

import pandas as pd
import numpy as np
import collections
import itertools
import random
import math
import sys
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


def peripheral(cells, **kwargs):
    """
    Peripheral metric
    """
    #handle kwargs
    processes = kwargs.get('processes',2)

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


def radial(cells, **kwargs):
    """
    Radial metric
    """
    #NOTE just returning the test metric for now
    return _test(cells, **kwargs)


def punctate(cells, **kwargs):
    """
    Punctate metric

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

    #NOTE make this parallel
    for cell_num,cell in enumerate(cells):
        if cell_num%50 == 0:
            sys.stdout.write('Punctate scoring cell {}\n'.format(cell_num))
            sys.stdout.flush()

        cell = cell.filter_genes_by_count(min_gene_spots=2)

        all_genes = np.array([g for z in cell.zslices for g in cell.spot_genes[z]])
        all_spots = np.array([xy for z in cell.zslices for xy in cell.spot_coords[z]])

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
            null_var = np.var([(exp-d)/exp for d in null_dists])

            data['metric'].append('puncta')
            data['cell_id'].append(cell.cell_id)
            data['annotation'].append(cell.annotation)
            data['num_spots'].append(cell.n)
            data['gene'].append(gene)
            data['num_gene_spots'].append(cell.gene_counts[gene])
            data['score'].append(score)
            data['variance'].append(null_var)


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

            data['metric'].append('periph')
            data['cell_id'].append(cell.cell_id)
            data['annotation'].append(cell.annotation)
            data['num_spots'].append(cell.n)
            data['gene'].append(gene)
            data['num_gene_spots'].append(cell.gene_counts[gene])
            data['score'].append(score)
            data['variance'].append(cell.gene_vars[gene])

    return pd.DataFrame(data)





def _peripheral_dist_and_rank(cell):
    """
    Helper function to calculate peripheral ranks
    """
    pass

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



def _central_dist_and_rank(cell):
    """
    Helper function to calculate central ranks
    """
    pass

