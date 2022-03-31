import SRRS
from SRRS import scoring

import sys

def main():

    hdf5_path = sys.argv[1]
    out_path = hdf5_path+'.shrunk_0.8'

    sample = SRRS.HDF5(hdf5_path)

    ######################
    #                    #
    #  shrink the cells  #
    #                    #
    ######################
    shrunk_cells = (c.shrink_boundaries(scale_factor=0.8) for c in sample.iter_cells())

    ################################################################################
    #                                                                              #
    #   re-filter since a cell might now not have enough RNA spots to be included  #
    #                                                                              #
    ################################################################################
    #Filter out cells with fewer than 100 total spots
    min_tot_spots = 100
    shrunk_cells = (c for c in shrunk_cells if c.n >= min_tot_spots)

    #Filter out cells with fewer than 10 unique genes that each have at least 3 spots
    min_unique_genes = 10
    min_spots_per_gene = 3
    shrunk_cells = (c for c in shrunk_cells
        if sum(v >= min_spots_per_gene for v in c.gene_counts.values()) >= min_unique_genes
    )


    ###############################################
    #                                             #
    #  recalculate the gene vars after shrinking  #
    #                                             #
    ###############################################
    shrunk_cells = scoring._iter_vars(shrunk_cells, processes=10)
    SRRS.HDF5.write_cells(shrunk_cells, out_path)


if __name__ == '__main__':
    main()

