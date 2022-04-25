import shapely.geometry

import pandas as pd
import numpy as np
import collections
import random
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


def punctate(cell, num_iterations=100):
    """
    Punctate metric

    Performs a first iterative step of randomly choosing two spots from each gene
    Then ranks the genes based on the distance between their respective two spots?
    Doesn't make sense to rank individual spots? since spots of the same gene will get the same rank?

    gene/cell score is the normalized average gene rank over all iterations
    Non-normalized score is expected to be (num_genes+1)/2
    """
    spot_gene_coords = [
        (g,(x,y)) for z in cell.zslices 
        for g,(x,y) in zip(cell.spot_genes[z],cell.spot_coords[z])
    ]

    all_ranks_df = pd.DataFrame()

    for i in range(num_iterations):
        #shuffle the spot order
        random.shuffle(spot_gene_coords)

        #choose two spots from each gene to be representatives (NOTE can do this more efficiently!)
        g_reps = {g:[] for g in cell.genes}
        for gene,xy in spot_gene_coords:
            if len(g_reps[gene]) < 2:
                g_reps[gene].append(xy)

        dists = []
        for gene in g_reps:
            x1,y1 = g_reps[gene][0]
            x2,y2 = g_reps[gene][1]
            dists.append((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
            
        gene_ranks = np.array(dists).argsort().argsort()+1 #add one so ranks start at 1 rather than 0

        #dists and ranks will always be in the order of the cell.genes
        gene_ranks_df = pd.DataFrame({'gene':cell.genes,'rank':gene_ranks,'iteration':i})
    
        all_ranks_df = pd.concat((all_ranks_df, gene_ranks_df),ignore_index=True)

    all_ranks_df['spot_count'] = all_ranks_df['gene'].map(cell.gene_counts)

    return all_ranks_df


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



