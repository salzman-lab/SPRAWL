import pandas as pd
import numpy as np
import collections
import random

from . import scoring
from . import metrics
from . import utils

def gene_celltype_sim_null(cells, metric, within_z=True, n_its=1000):
    """
    Perform a null simulation on cells by randomly changing gene labels

    For each gene/celltype pair, calculate how many have more extreme scores than
    expected by chance

    Steps
    0. For each iteration
    1. Permute gene labels
    2. Calculate SRRS scores per gene/cell using metric
    3. Group cells by cell-type
    4. Count how many permutations resulted in significant peripheral distributions

    yield pandas dataframes of gene_celltype_scoring over the multiple iterations
    """

    cells = list(utils._iter_vars(cells)) #avoid consuming iterator

    #iterate through the permuatations
    for it_num in range(n_its):

        #permute all the cell gene labels (occurs in-place)
        for cell in cells:
            null_permute_gene_labels(cell, within_z)

        #score each gene/cell pair
        srrs_df = pd.concat(scoring.iter_scores(cells, metric))

        #calculate stats for each gene/celltype pair
        agg_df = scoring.gene_celltype_scoring(srrs_df)
        agg_df['it_num'] = it_num

        yield agg_df



def gene_cell_sim_null(cells, metric, within_z=True, n_its=1000, alpha=0.05, processes=2):
    """
    Perform a null simulation on cells by randomly changing gene labels

    For each gene/cell pair, calculate how many permutations have more extreme scores than
    expected by chance

    Return pandas dataframe with the following columns
    * cell_id
    * gene
    * perms_lt_score
    * perms_gt_score
    * sig
    """

    obs_score_df = scoring.iter_scores(cells, metric=metric, processes=processes)
    obs_scores = obs_score_df.set_index(['cell_id','gene'])['score'].sort_index()
    counts_lt = pd.Series(data=0.0, index=obs_scores.index, name='perms_lt_score')

    #for each iteration, permute gene labels and count how many permutations have lower scores than obs
    for _ in range(n_its):
        for cell in cells:
            null_permute_gene_labels(cell, within_z)

        perm_score_df = scoring.iter_scores(cells, metric=metric, processes=processes)
        perm_scores = perm_score_df.set_index(['cell_id','gene'])['score'].sort_index()

        counts_lt += perm_scores < obs_scores
       
    #calculate two-sided probabilities
    sig_cutoff = 0.5*n_its*alpha
    df = counts_lt.reset_index()
    df['perms_gt_score'] = n_its-df['perms_lt_score']
    df['sig'] = df['perms_lt_score'].le(sig_cutoff) | df['perms_gt_score'].le(sig_cutoff)
 
    return df


def null_permute_gene_labels(cell, within_z=True):
    """
    Take as input a Cell object
    permutes in place but also returns Cell object reference

    doesn't matter which metric is used

    within_z=True means gene labels will only be reassigned within each z-slice
    within_z=False means gene labels will be reassigned cell-wide
    """
    if within_z:
        for z in cell.zslices:
            np.random.shuffle(cell.spot_genes[z])

    else:
        all_genes = [g for z in cell.zslices for g in cell.spot_genes[z]]
        random.shuffle(all_genes)

        i = 0
        for z in cell.zslices:
            slice_spots = len(cell.spot_genes[z])
            cell.spot_genes[z] = all_genes[i:i+slice_spots]
            i += slice_spots

    return cell



