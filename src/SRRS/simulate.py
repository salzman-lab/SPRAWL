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

    cells = list(scoring._iter_vars(cells)) #avoid consuming iterator

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
    obs_scores = obs_score_df.set_index(['cell_id','gene'])['score']
    counts_lt = pd.Series(data=0, index=obs_scores.index, name='perms_lt_score')

    #for each iteration, permute gene labels and count how many permutations have lower scores than obs
    for _ in range(n_its):
        for cell in cells:
            null_permute_gene_labels(cell, within_z)

        perm_score_df = scoring.iter_scores(cells, metric=metric, processes=processes)
        perm_scores = perm_score_df.set_index(['cell_id','gene'])['score']
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


def lin_permute_gene_labels(cell, metric, gene_ks={}):
    """
    Take as input a Cell object
    permutes in place but also returns Cell object reference

    Assigns gene labels biased by spot ranks and gene_ks dict
    under the assumption that the probability distribution is linear

    k represents the ratio in of p(g = rank1)/p(g = rankN)

    if a gene is not present in a gene_ks, then a uniform distribution is assumed (k = 1)
    """

    #calculate spot ranks
    metric_f = scoring.available_metrics[metric]
    cell = metric_f(cell)
    n = cell.n

    #keep track of how many spots of each gene we need to assign
    remaining_gene_counts = cell.gene_counts.copy()

    #create a lookup of normalized probabilities for each rank
    p_gene_ranks = {}
    for gene in remaining_gene_counts.keys():
        k = gene_ks[gene] if gene in gene_ks else 1

        if k != 1:
            step = (1-k)/(n-1)
            ps = np.arange(k,1+step,step)
            p_gene_ranks[gene] = ps/sum(ps)
        else:
            p_gene_ranks[gene] = [1/n]*n

    #Update the median ranks after permuting the gene labels
    cell = _alt_perm_helper(cell, p_gene_ranks)
    metrics._update_med_ranks(cell)
    return cell


def exp_permute_gene_labels(cell, metric, gene_ks={}):
    """
    Take as input a Cell object
    permutes in place but also returns Cell object reference

    Assigns gene labels biased by spot ranks and gene_ks dict
    under the assumption that the probability distribution follows an exponential
    function with equation: p = exp(-k*rank/n)

    large magnitude positive values of k are highly biased towards smaller ranks
    large magnitude negative values of k are highly biased towards larger ranks

    if a gene is not present in a gene_ks, then a uniform distribution is assumed (k = 0)
    """
    #calculate spot ranks
    metric_f = scoring.available_metrics[metric]
    cell = metric_f(cell)
    n = cell.n

    #create a lookup of normalized probabilities for each rank
    p_gene_ranks = {}
    for gene in cell.gene_counts.keys():
        k = gene_ks[gene] if gene in gene_ks else 0
        ranks = np.arange(1,n+1)
        ps = np.exp(-k*ranks/n)
        p_gene_ranks[gene] = ps/sum(ps)

    #Update the median ranks after permuting the gene labels
    cell = _alt_perm_helper(cell, p_gene_ranks)
    metrics._update_med_ranks(cell)
    return cell


def _alt_perm_helper(cell, p_gene_ranks):
    """
    Helper to iterate through all z's and all RNA spots to re-assign genes
    """
    #keep track of how many spots of each gene we need to assign
    remaining_gene_counts = cell.gene_counts.copy()

    for z in cell.zslices:
        new_spot_genes = []
        for spot_rank in cell.spot_ranks[z]:
            #calculate normalized p's weights of choosing each gene
            genes = []
            raw_weights = []
            for gene,count in remaining_gene_counts.items():
                p_gene = p_gene_ranks[gene][spot_rank-1] #ranks are 1-indexed
                genes.append(gene)
                raw_weights.append(p_gene*count)

            #sometimes all the probabilities will be too small and rounded to 0
            #in this case, choose uniformly between the genes depending only on remaining count
            if sum(raw_weights) == 0:
                raw_weights = remaining_gene_counts.values()

            weights = [w/sum(raw_weights) for w in raw_weights]

            #choose a gene based on the weights
            gene = np.random.choice(genes,p=weights)

            #Add chosen gene to list
            new_spot_genes.append(gene)

            #update remaining gene counts
            if remaining_gene_counts[gene] > 1:
                remaining_gene_counts[gene] -= 1
            else:
                del remaining_gene_counts[gene]

        #Assign chosen genes
        cell.spot_genes[z] = new_spot_genes

    #Update the median ranks after permuting the gene labels
    metrics._update_med_ranks(cell)
    return cell


