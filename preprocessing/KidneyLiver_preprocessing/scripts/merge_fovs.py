#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import h5py
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Merge HDF5s together')

    parser.add_argument('--hdf5_fovs', dest='hdf5_fovs', required=True, nargs='+',
                        help='Paths to the HDF5 different fovs to merge')

    parser.add_argument('--sample', dest='sample', required=True,
                        help='Just used to name combined output')

    args = parser.parse_args()
    return args


# Main #
if __name__ == '__main__':
    args = parse_args()

    #Create the combined HDF5 for final output
    out_f = h5py.File('{}.hdf5'.format(args.sample),'w')
    out_cells_grp = out_f.create_group('cells')

    #Iterate through the fovs adding all cells to a new output hdf5
    unique_cell_ids = set()
    unique_genes = set()

    for in_hdf5_path in args.hdf5_fovs:
        in_hdf5 = h5py.File(in_hdf5_path,'r')

        for cell_id in in_hdf5['cells']:
            if cell_id in unique_cell_ids:
                sys.stderr.write('ERROR found cell {} again!\n'.format(cell_id))
                sys.exit(1)

            in_hdf5.copy('cells/{}'.format(cell_id), out_cells_grp)

            #keep track of the cell_ids
            unique_cell_ids.add(cell_id.encode())

        unique_genes = unique_genes.union(in_hdf5['genes'])

    out_f.create_dataset('cell_ids',data=sorted(list(unique_cell_ids)))
    out_f.create_dataset('genes',data=sorted(list(unique_genes)))

    #Close the output file
    out_f.close()

