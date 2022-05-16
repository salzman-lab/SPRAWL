import SRRS

import pandas as pd
import sys
import os

def main():
    hdf5_path = '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/inputs/mouse1sample1.hdf5'
    out_name = 'refactored_mouse1sample1.csv'

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

    #score and write out
    scores_df = SRRS.iter_scores(cells, metric='peripheral')
    scores_df.to_csv(out_name, index=False)

if __name__ == '__main__':
    main()

