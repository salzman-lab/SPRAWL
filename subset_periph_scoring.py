import SRRS
from SRRS import scoring

import pandas as pd
import time
import sys
import os

def main():
    prog_start = time.time()
    hdf5_path = 'inputs/m1s1_subset.hdf5'

    sample = SRRS.HDF5(hdf5_path)
    cells = sample.iter_cells()

    #cache the expensive var calculations without filtering any cells out
    cells = list(scoring._iter_vars(cells))
    sample.save_gene_vars(cells)


if __name__ == '__main__':
    main()

