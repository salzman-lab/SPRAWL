import SRRS

import pandas as pd
import time
import sys
import os

def main():
    prog_start = time.time()
    hdf5_path = sys.argv[1]
    out_name = os.path.basename(hdf5_path).split('.')[0]+'_no_drop_periph_scores.csv'

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
    scores = SRRS.iter_scores(cells, metric='peripheral')

    for i,score_df in enumerate(scores):
        score_df.to_csv(out_name, mode='a', index=False, header = i==0)

if __name__ == '__main__':
    main()

