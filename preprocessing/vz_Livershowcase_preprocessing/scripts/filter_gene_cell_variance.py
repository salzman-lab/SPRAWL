#!/usr/bin/env python3
import sys
sys.path.append('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/src') #ugly

import SRRS
from SRRS import scoring

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Filter cells and calculate gene/cell variances')

    parser.add_argument('--hdf5_path', dest='hdf5_path', required=True,
                        help='Path to the hdf5 file to be processed')

    parser.add_argument('--out_path', dest='out_path', required=True,
                        help='Path to output hdf5 file')

    parser.add_argument('--min_gene_spots', dest='min_gene_spots', required=True, type=int,
                        help='Minimum spots for a gene to be kept in a cell')

    parser.add_argument('--min_tot_spots', dest='min_tot_spots', required=True, type=int,
                        help='Minimum spots in a cell to be filtered-in')

    parser.add_argument('--min_genes', dest='min_genes', required=True, type=int,
                        help='Minimum unique genes in a cell to be filtered-in')


    args = parser.parse_args()
    return args




def main():
    args = parse_args()

    sample = SRRS.HDF5(args.hdf5_path)
    cells = sample.iter_cells()

    #Filter out genes in cells that have only 1 spot
    cells = (c.filter_low_count_genes(args.min_gene_spots) for c in cells)

    #Then filter out cells with fewer than threshold total spots
    cells = (c for c in cells if c.n >= args.min_tot_spots)

    #Filter out cells with fewer than threshold unique genes
    cells = (c for c in cells if len(c.genes) >= args.min_genes)

    #cache the expensive var calculations and write out to new hdf5
    cells = scoring._iter_vars(cells, processes=10) #NOTE
    SRRS.HDF5.write_cells(cells, args.out_path)


if __name__ == '__main__':
    main()

