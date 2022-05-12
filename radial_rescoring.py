import SRRS
from SRRS import simulate

import pandas as pd
import sys
import os

def main():
    hdf5_path = sys.argv[1]

    metric = 'radial'

    out_name = os.path.basename(hdf5_path).split('.')[0]+'_{}_permuted_scores.csv'.format(metric)

    sample = SRRS.HDF5(hdf5_path)
    cells = sample.iter_cells()

    #Filter out cells with fewer than 100 total spots
    min_tot_spots = 100
    cells = (c for c in cells if c.n >= min_tot_spots)

    #Filter out cells with fewer than 10 unique genes that each have at least 3 spots
    min_unique_genes = 10
    min_spots_per_gene = 3
    cells = (c for c in cells
        if sum(v >= min_spots_per_gene for v in c.gene_counts.values()) >= min_unique_genes
    )

    #Permute the gene labels in all cells to make a null dist
    cells = (simulate.null_permute_gene_labels(cell, within_z=False) for cell in cells)

    #score and write out
    scores_df = SRRS.iter_scores(cells, metric=metric, num_pairs=10, processes=10)
    scores_df.to_csv(out_name, index=False)

if __name__ == '__main__':
    main()

