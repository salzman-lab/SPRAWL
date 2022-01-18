#!/usr/bin/env python3
from rtree import index
import pandas as pd
import numpy as np
import argparse

import shapely.geometry
import h5py

def parse_args():
    parser = argparse.ArgumentParser(description='Assign spots to cells')

    parser.add_argument('--segm_path', dest='segm_path', required=True,
                        help='Path to the segment file to be processed')

    parser.add_argument('--rtree_name', dest='rtree_name', required=True,
                        help='Name of the RTree of the spots')

    parser.add_argument('--z', dest='z', required=True, type=int,
                        help='Which z-slice to process')

    args = parser.parse_args()
    return args


def make_polygons(r):
    xs = [float(x) for x in r[0].split(', ')]
    ys = [float(y) for y in r[1].split(', ')]
    
    #Should be the same number of xs and ys
    if len(xs) != len(ys):
        print('Different number of coords for x,y')
        return None

    #Need a minimum of 3 points to make a polygon
    if len(xs) < 3:
        return None
    
        
    poly = shapely.geometry.Polygon(zip(xs,ys))
    return poly
    

# Main #
if __name__ == '__main__':
    args = parse_args()

    #Convert the cells into polygons
    seg_df = pd.read_csv(args.segm_path, index_col=0)

    x_col = 'boundaryX_z{}'.format(args.z)
    y_col = 'boundaryY_z{}'.format(args.z)

    cell_polys = seg_df[[x_col, y_col]].dropna().apply(make_polygons,axis=1).dropna()

    #Read in the RTree of the spots
    idx = index.Rtree(args.rtree_name)

    #Create the HDF5 file to prepare to write out
    out_f = h5py.File('slice_cells.hdf5','w')
    out_f.attrs['zslice'] = args.z

    cells_grp = out_f.create_group('cells')
    all_genes = set()
    cell_ids = []

    #Loop through each cell 
    for cell_id,poly in cell_polys.iteritems():

        cell_genes = []
        cell_spots = []

        #Find which spots are in this cell
        boundary_hits = idx.intersection(poly.bounds, objects=True)
        
        for hit in boundary_hits:
            gene = hit.object
            x,y = hit.bbox[:2]
        
            point = shapely.geometry.Point(x,y)
            
            if poly.contains(point):        
                cell_genes.append(gene.encode())
                cell_spots.append([x,y])
      

        #Skip cell-slice if no RNA spots are present
        if len(cell_spots) == 0:
            continue
 
        #Instantiate hdf5 groups
        cell_grp = cells_grp.create_group(cell_id)

        cell_grp.create_dataset('boundary', data=poly.exterior.coords)
        cell_grp.create_dataset('spot_genes', data=cell_genes)
        cell_grp.create_dataset('spot_coords', data=cell_spots)

        #Update the overall cell-ids list and the overall genes
        cell_ids.append(cell_id.encode())
        all_genes = all_genes.union(cell_genes)


    #Create datasets for the cell_ids and genes at the root level
    out_f.create_dataset('cell_ids',data=cell_ids)
    out_f.create_dataset('genes',data=sorted(list(all_genes)))


    #Close the HDF5 file
    out_f.close()


