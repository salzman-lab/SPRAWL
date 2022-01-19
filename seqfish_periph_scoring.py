import SRRS

import pandas as pd
import time
import sys
import os

def main():
    prog_start = time.time()
    hdf5_path = 'preprocessing/SeqFishplus_preprocessing/seqfish_plus.hdf5'
    out_name = 'seqfish_peripheral_gene_cell.csv'

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

    #cache the expensive var calculations
    cells = list(scoring._iter_vars(cells))
    sample.save_gene_vars(cells)

    #score and write out
    scores = SRRS.iter_scores(cells, metric='peripheral')
    for i,score_df in enumerate(scores):
        score_df.to_csv(out_name, mode='a', index=False, header = i==0)


if __name__ == '__main__':
    main()

