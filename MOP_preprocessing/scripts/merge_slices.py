#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
import h5py

def parse_args():
    parser = argparse.ArgumentParser(description='Merge HDF5s together')

    parser.add_argument('--slices', dest='slices', required=True, nargs='+',
                        help='Paths to the HDF5 z-slices')

    parser.add_argument('--cell_labels', dest='cell_labels', required=True,
                        help='Path to the cell-type annotation table')

    parser.add_argument('--mouse_sample', dest='mouse_sample', required=True,
                        help='Just used to name combined output')

    args = parser.parse_args()
    return args


# Main #
if __name__ == '__main__':
    args = parse_args()
    print(args)


    #Create the combined HDF5 for final output
    out_f = h5py.File('{}.hdf5'.format(args.mouse_sample),'w')


    #Make a dict of cell_id --> annotation
    ann_df = pd.read_csv(
        args.cell_labels,
        index_col=0,
    )
    cell_id_to_ann = dict(ann_df['label'])

    #Read in the hdf5s per slice
    slice_hdf5s = [h5py.File(sl,'r') for sl in args.slices]

    unique_cell_ids = list(set(cell_id for f in slice_hdf5s for cell_id in f['cell_ids']))
    unique_genes = sorted(list(set(gene for f in slice_hdf5s for gene in f['genes'])))

    out_f.create_dataset('cell_ids',data=unique_cell_ids)
    out_f.create_dataset('genes',data=unique_genes)

    cells_grp = out_f.create_group('cells')

    for cell_id in unique_cell_ids:
        cell_grp = cells_grp.create_group(cell_id)
        
        boundaries = cell_grp.create_group('boundaries')
        spot_coords = cell_grp.create_group('spot_coords')
        spot_genes = cell_grp.create_group('spot_genes')
        
        zs = []
        cell_genes = set()
        num_cell_spots = 0
        
        for slice_hdf5 in slice_hdf5s:
            
            if cell_id not in slice_hdf5['cells']:
                continue
                
            z = str(slice_hdf5.attrs['zslice'])
            zs.append(z)
            
            sl_cell = slice_hdf5['cells'][cell_id]
            sl_spots = sl_cell['spot_coords'][:]
            sl_genes = sl_cell['spot_genes'][:]
            
            boundaries.create_dataset(z, data=sl_cell['boundary'][:])
            spot_coords.create_dataset(z, data=sl_spots)
            spot_genes.create_dataset(z, data=sl_genes)
            
            num_cell_spots += len(sl_spots)
            cell_genes = cell_genes.union(sl_genes)
        
        cell_grp.attrs['zslices'] = sorted(zs)
        cell_grp.attrs['num_genes'] = len(cell_genes)
        cell_grp.attrs['num_spots'] = num_cell_spots
        cell_grp.attrs['annotation'] = cell_id_to_ann.get(cell_id.decode(), 'unannotated')
        
    for f in slice_hdf5s:
        f.close()


