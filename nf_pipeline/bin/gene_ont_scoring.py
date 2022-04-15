#!/usr/local/bin/python
#Above is python path in the docker container with SRRS installed
import SRRS
from SRRS import scoring

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene_cell_path')
    parser.add_argument('--output_name')
    parser.add_argument('--min_cells_per_gene_ont', type=int)
    args = parser.parse_args()

    agg_df = scoring.gene_celltype_scoring(
        args.gene_cell_path,
        min_cells_per_gene_ont=args.min_cells_per_gene_ont,
        extra_cols = {'':''},
    )

    agg_df.to_csv(args.output_name, index=False)

    

if __name__ == '__main__':
    main()

