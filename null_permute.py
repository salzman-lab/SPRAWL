import SRRS
from SRRS import scoring,simulate

import sys

def main():
    hdf5_path = 'inputs/mouse1sample3.hdf5'
    out_name = 'm1s3_radial_null_permutes_across_z.csv'

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

    #Filter to maintain only Pvalb and VIP celltypes and subtypes
    cells = (c for c in cells
        if ('Pvalb' in c.annotation) or ('Vip' in c.annotation)
    )

    #calculate and save out csv of simulations
    sim_df_iter = simulate.gene_celltype_sim_null(cells, 'radial', within_z=False, n_its=1000)

    for i,sim_df in enumerate(sim_df_iter):
        sim_df.to_csv(out_name, mode='a', index=False, header = i==0)
    


if __name__ == '__main__':
    main()

