#!/usr/bin/env python3
from rtree import index
import pandas as pd
import numpy as np
import argparse
import collections

import shapely.geometry
import h5py

def parse_args():
    parser = argparse.ArgumentParser(description='Assign spots to cells')

    parser.add_argument('--hdf5_path', dest='hdf5_path', required=True,
                        help='Path to the hdf5 file to be processed')

    parser.add_argument('--rtree_name', dest='rtree_name', required=True,
                        help='RTree name')

    args = parser.parse_args()
    return args



# Main #
if __name__ == '__main__':
    args = parse_args()
    in_f = h5py.File(args.hdf_path,'r')

    #Create the HDF5 file to prepare to write out
    out_f = h5py.File('spot_assigned_fov.hdf5','w')
    out_cells_grp = out_f.create_group('cells')
    all_cell_ids = []
    all_genes = set()

    #Read in the 3D RTree of the spots
    p = index.Property()
    p.dimension = 3
    idx = index.Rtree(args.rtree_name, properties=p)

    for cell_id in in_f['featuredata']:
        in_cell = in_f['featuredata'][cell_id]
        min_x,min_y,max_x,max_y = in_cell.attrs['bounding_box']

        cell_genes = collections.defaultdict(list)
        cell_spots = collections.defaultdict(list)

        polys = {}

        for z,_ in enumerate(in_cell['z_coordinates']):

            #Find which spots are in this cell boundary box
            boundary_hits = idx.intersection((min_x, min_y, z, max_x, max_y, z), objects=True)

            #Create a Shapely polygon from the x,y coords at the given z-slice
            z_slice = in_cell['zIndex_{}'.format(z)]
            if not len(cell[z_slice]):
                continue

            xy_s = in_cell[z_slice]['p_0']['coordinates'][0,:,:]
            poly = shapely.geometry.Polygon()
            polys[z] = poly

            for hit in boundary_hits:
                gene = hit.object
                x,y = hit.bbox[:2]
            
                point = shapely.geometry.Point(x,y)
                
                if poly.contains(point):        
                    cell_genes[z].append(gene.encode()) #TODO need to convert barcode_id to gene name
                    cell_spots[z].append([x,y])
          
        #Skip creating cell if no RNA spots are present in any slice
        if len(cell_spots) == 0:
            continue

        #Instantiate hdf5 group for the cell
        cell_grp = out_cells_grp.create_group(cell_id)
        spot_genes_grp = cell_grp.create_group('spot_genes')
        spot_coords_grp = cell_grp.create_group('spot_coords')
        boundary_grp = cell_grp.create_group('boundaries')

        uniq_cell_genes = {g for z,gs in cell_genes.items() for g in gs})

        #Add attrs to the cell
        cell_grp.attrs['annotation'] = 'None' #TODO add annotation data to the cells
        cell_grp.attrs['num_spots'] = sum(len(s) for z,s in cell_spots.items())
        cell_grp.attrs['num_genes'] = len(uniq_cell_genes)
        cell_grp.attrs['zslices'] = [str(z) for z in cell_spots.keys()]

        #Loop through the z-slices adding the spot-gene, spot-coords, and boundaries
        for z,coords in cell_spots.items():
            spot_coords_grp.create_dataset(z, data=coords)

        for z,genes in cell_genes.items():
            spot_genes_grp.create_dataset(z, data=genes)

        for z,poly in polys.items():
            boundary_grp.create_dataset(z, data=poly.exterior.coords)

        #Update the overall cell-ids list and the overall genes
        all_cell_ids.append(cell_id.encode())
        all_genes = all_genes.union(uniq_cell_genes)


    #Create datasets for the cell_ids and genes at the root level
    out_f.create_dataset('cell_ids',data=all_cell_ids)
    out_f.create_dataset('genes',data=sorted(list(all_genes)))

    #Close the HDF5 files
    in_f.close()
    out_f.close()


