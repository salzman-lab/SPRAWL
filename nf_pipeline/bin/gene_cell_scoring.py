#!/usr/local/bin/python
#Above is python path in the docker container with SRRS installed
import SRRS
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hdf5_path')
    parser.add_argument('--metric')
    parser.add_argument('--experiment')
    parser.add_argument('--sample')
    parser.add_argument('--output_name')
    parser.add_argument('--min_genes_per_cell', type=int)
    parser.add_argument('--min_tot_counts_per_cell', type=int)
    args = parser.parse_args()

    sample = SRRS.HDF5(args.hdf5_path)
    cells = sample.iter_cells()

    #Filter out cells with fewer than 100 total spots
    cells = (
        c for c in cells
        if c.n >= args.min_tot_counts_per_cell and len(c.genes) >= args.min_genes_per_cell
    )

    #score and write out
    scores = SRRS.iter_scores(cells, metric=args.metric)

    for i,score_df in enumerate(scores):
        score_df['experiment'] = args.experiment
        score_df['sample'] = args.sample
        score_df.to_csv(args.output_name, mode=('w' if i == 0 else 'a'), index=False, header = i==0)

if __name__ == '__main__':
    main()

