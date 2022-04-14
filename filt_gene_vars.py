import SRRS
from SRRS import scoring

import pandas as pd
import time
import sys
import os

def main():
    prog_start = time.time()

    hdf5_path = sys.argv[1]
    out_path = hdf5_path+'.gv_filt.gene_count_filt'

    sample = SRRS.HDF5(hdf5_path)
    cells = sample.iter_cells()

    #Filter out genes in cells that have only 1 spot
    min_gene_spots = 2
    cells = (c.filter_low_count_genes(min_gene_spots) for c in cells)

    #Then filter out cells with fewer than 100 total spots
    min_tot_spots = 100
    cells = (c for c in cells if c.n >= min_tot_spots)

    #Filter out cells with fewer than 10 unique genes that each have at least 3 spots
    min_unique_genes = 10
    min_spots_per_gene = 3
    cells = (c for c in cells
        if sum(v >= min_spots_per_gene for v in c.gene_counts.values()) >= min_unique_genes
    )

    #cache the expensive var calculations and write out to new hdf5
    cells = scoring._iter_vars(cells)
    SRRS.HDF5.write_cells(cells, out_path)


if __name__ == '__main__':
    main()

