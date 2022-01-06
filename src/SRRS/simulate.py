import pandas as pd
import numpy as np
import collections
import random

from . import scoring
from . import metrics
from . import utils

def gene_celltype_sim_null(cells, metric, within_z=True, n_its=1000, alpha=0.05):
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

    Return pandas dataframe with the following columns
    * metric
    * gene
    * annotation
    * median_gene_spots
    * median_total_spots
    * num_cells
    * alpha
    * z_score
    * two_sided_p

    """

    cells = list(scoring._iter_vars(cells)) #avoid consuming iterator

    data = {
        'metric':[],
        'gene':[],
        'annotation':[],
        'median_gene_spots':[],
        'median_total_spots':[],
        'num_cells':[],
        'alpha':[],
        'z_score':[],
        'two_sided_p':[],
    }

    #iterate through the permuatations
    for _ in range(n_its):

        #permute all the cell gene labels
        for cell in cells:
            null_permute_gene_labels(cell, within_z)

        #score
        scoring_df = pd.concat(scoring.iter_scores(cells, metric))




def gene_cell_sim_null_peripheral(cells, within_z=True, n_its=1000, alpha=0.05):
    """
    Perform a null simulation on cells by randomly changing gene labels

    For each gene/cell pair, calculate how many have more extreme scores than
    expected by chance

    Steps
    1. Permute gene labels
    2. Calculate SRRS scores using metric
    3. Count how many genes have significant peripheral distributions

    Return pandas dataframe with the following columns
    * metric
    * cell_id
    * annotation
    * num_spots
    * gene
    * num_gene_spots
    * variance
    * num_its
    * alpha
    * num_sig_its

    Calculate theoretical variance just once
    Calculate spot ranks just once, then permute spots
    """
    meds_per_cell_per_gene = {}
    cells = list(scoring._iter_vars(cells)) #avoid consuming generator

    #make all the permutations and store the results
    m_n_meds = collections.defaultdict(set)
    for cell in cells:
        meds_per_cell_per_gene[cell.cell_id] = collections.defaultdict(list)
        n = cell.n
        _,ranks = metrics._peripheral_dist_and_rank(cell)
        for _ in range(n_its):
            null_permute_gene_labels(cell, within_z)
            genes = np.array([g for z in cell.zslices for g in cell.spot_genes[z]])

            for g,m in cell.gene_counts.items():
                med = np.median(ranks[genes == g])+1
                meds_per_cell_per_gene[cell.cell_id][g].append(med)
                m_n_meds[(m,n)].add(med)

    #multiprocess the p_twosided calcs which are slow
    p_cache = scoring._calc_p_twosided(m_n_meds)


    #consolidate data
    data = {
        'metric':'peripheral',
        'cell_id':[],
        'annotation':[],
        'num_spots':[],
        'gene':[],
        'num_gene_spots':[],
        'theory_variance':[],
        'emp_variance':[],
        'num_its':[],
        'alpha':[],
        'num_sig_its':[],
    }

    for cell in cells:
        for g,m in cell.gene_counts.items():
            obs_meds = meds_per_cell_per_gene[cell.cell_id][g]

            num_sig = sum(p_cache[(m,cell.n,med)] < alpha for med in obs_meds)

            data['cell_id'].append(cell.cell_id)
            data['annotation'].append(cell.annotation)
            data['num_spots'].append(cell.n)
            data['gene'].append(g.decode())
            data['num_gene_spots'].append(m)
            data['theory_variance'].append(cell.gene_vars[g])
            data['emp_variance'].append(np.var(obs_meds))
            data['num_its'].append(n_its)
            data['alpha'].append(alpha)
            data['num_sig_its'].append(num_sig)

    return pd.DataFrame(data)


def null_permute_gene_labels(cell, within_z=True):
    """
    Take as input a Cell object
    permutes in place but also returns Cell object reference

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


