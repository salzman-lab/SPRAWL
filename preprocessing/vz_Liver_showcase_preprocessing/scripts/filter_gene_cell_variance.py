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

    parser.add_argument('--min_spots', dest='min_spots', required=True, type=int,
                        help='Minimum spots in a cell to be filtered-in')

    parser.add_argument('--min_genes', dest='min_genes', required=True, type=int,
                        help='Minimum unique genes in a cell to be filtered-in')


    args = parser.parse_args()
    return args




def main():
    args = parse_args()

    sample = SRRS.HDF5(args.hdf5_path)

    #Filter out cells with fewer than X total spots
    #Filter out cells with fewer than Y unique genes
    cells = (
        c for c in sample.iter_cells()
        if c.n >= args.min_spots and len(c.genes) >= args.min_genes
    )

    #cache the expensive var calculations and write out to new hdf5
    cells = scoring._iter_vars(cells, processes=10) #NOTE
    SRRS.HDF5.write_cells(cells, args.out_path)


if __name__ == '__main__':
    main()

