import pandas as pd
import numpy as np
import random

from . import scoring

def sim_null(cells, metric, n_its=100):
    """
    Perform a null simulation on cells by randomly changing gene labels

    Steps
    1. Permute gene labels
    2. Calculate SRRS scores using metric
    3. Count how many genes are

    Return pandas dataframe with the following columns
    * metric
    * cell_id
    * annotation
    * num_spots
    * gene
    * num_gene_spots
    * variance
    * num_sig_its

    Avoid re-calculating variance since its slow and unnecessary
    """
    data = {
        'metric':metric,
        'cell_id':[],
        'annotation':[],
        'num_spots':[],
        'gene':[],
        'num_gene_spots':[],
        'variance':[],
        'num_sig_its':[],
    }

    cells = scoring._iter_vars(cells)

    for cell in cells:
        data['cell_id'].append(cell.cell_id)
        data['annotation'].append(cell.annotation)
        data['num_spots'].append(cell.n)

        sig_its_by_gene = {g:0 for g in cell.genes}

        for _ in range(n_its):
            null_permute_gene_labels(cell)
            #score

            for g,med in cell.gene_med_ranks.items():
                p_med = 

        for g,c in cell.gene_counts.items():
            data['gene'].append(g.decode())
            data['num_gene_spots'].append(c)
            data['variance'].append(cell.gene_vars[g])
            data['num_sig_its'].append(sig_its_by_gene[g])

    return pd.DataFrame(data)


def null_permute_gene_labels(cell, within_z=True):
    """
    Take as input a Cell object
    no return, permutes in place

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


