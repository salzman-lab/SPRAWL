
import SRRS
from SRRS import simulate

import sys
import os

def main():
    hdf5_path = sys.argv[1]
    metric = sys.argv[2]

    out_stem = os.path.basename(hdf5_path).split('.')[0]+'_{}'.format(metric)

    #Read in the cells
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

    #Permute the gene labels of every cell
    cells = (simulate.null_permute_gene_labels(c, within_z=False) for c in cells)

    #gene/cell score
    scores_df = SRRS.iter_scores(cells, metric=metric)
    scores_df.to_csv(out_stem+'_gene_cells.csv', index=False)

    #gene/ont score
    gene_ont_df = SRRS.gene_celltype_scoring(scores_df, min_cells_per_gene_ont=10)

    #write gene/ont output 
    gene_ont_df.to_csv(out_stem+'_gene_onts.csv', index=False)


if __name__ == '__main__':
    main()

