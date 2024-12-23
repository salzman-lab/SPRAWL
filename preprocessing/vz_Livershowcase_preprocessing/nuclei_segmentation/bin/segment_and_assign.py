#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import h5py

from skimage import io
from skimage.segmentation import watershed, mark_boundaries
from skimage.filters import sobel
from skimage.exposure import histogram
from scipy import ndimage as ndi
from skimage.measure import find_contours, approximate_polygon
from shapely.geometry import Point, Polygon
import geopandas

import argparse
import sys
import os

#globally store entire image
im = None

def get_fov_boundaries(hdf5_fov):
    #Get the fov bounds in microns
    fov_micron_min_x = None
    fov_micron_min_y = None
    fov_micron_max_x = None
    fov_micron_max_y = None

    with h5py.File(hdf5_fov,'r') as f:
        for cell_id in f['featuredata']:
            micron_min_x,micron_min_y,micron_max_x,micron_max_y = f['featuredata'][cell_id].attrs['bounding_box']

            if not fov_micron_min_x or micron_min_x < fov_micron_min_x:
                fov_micron_min_x = micron_min_x
            if not fov_micron_min_y or micron_min_y < fov_micron_min_y:
                fov_micron_min_y = micron_min_y
            if not fov_micron_max_x or micron_max_x > fov_micron_max_x:
                fov_micron_max_x = micron_max_x
            if not fov_micron_max_y or micron_max_y > fov_micron_max_y:
                fov_micron_max_y = micron_max_y

    return fov_micron_min_x,fov_micron_min_y,fov_micron_max_x,fov_micron_max_y


def nuclei_segment_fov(transform, hdf5_fov):
    """
    Segment the nuclei within a field of view and assign nuclei to cells by boundary
    """
    #Read in the pixel to micron scaling and make helper functions
    A_micron_to_mosaic = np.genfromtxt(
        transform,
        delimiter=' ',
    )
    A_mosaic_to_micron = np.linalg.inv(A_micron_to_mosaic)#get the opposite direction transformation matrix

    pixel_to_micron = lambda x,y: np.matmul(A_mosaic_to_micron,[[x],[y],[1]])[:2].flatten()
    micron_to_pixel = lambda x,y: np.matmul(A_micron_to_mosaic,[[x],[y],[1]])[:2].flatten()

    #Get the boundaries to use to subset the image and convert them to pixels
    fov_micron_min_x,fov_micron_min_y,fov_micron_max_x,fov_micron_max_y = get_fov_boundaries(hdf5_fov)
    fov_pixel_min_x,fov_pixel_min_y = micron_to_pixel(fov_micron_min_x, fov_micron_min_y)
    fov_pixel_max_x,fov_pixel_max_y = micron_to_pixel(fov_micron_max_x, fov_micron_max_y)


    fov_img = im[
        int(fov_pixel_min_y):int(fov_pixel_max_y),
        int(fov_pixel_min_x):int(fov_pixel_max_x)
    ]

    #Use quantiles as the masking thresholds
    max_noise_thresh,min_signal_thresh = np.quantile(fov_img.flatten(),[0.2,0.95])

    elevation_map = sobel(fov_img)

    markers = np.zeros_like(fov_img)
    markers[fov_img < max_noise_thresh] = 1
    markers[fov_img > min_signal_thresh] = 2

    segmentation = watershed(elevation_map, markers)
    segmentation = ndi.binary_fill_holes(segmentation - 1)

    #THIS GIVES BACK (y,x) for some reason!!
    #followed tutorial here: https://scikit-image.org/docs/dev/auto_examples/edges/plot_polygon.html#sphx-glr-auto-examples-edges-plot-polygon-py
    raw_nuclei_polys = find_contours(segmentation,0)


    micron_nuclei_polys = []
    cell_ids = []

    #Iterate through returned nuclei polys and test if they fit within any of the cell boundaries
    with h5py.File(hdf5_fov,'r') as f:
        for cell_id in f['featuredata']:
            cell = f['featuredata'][cell_id]
            cell_micron_coords = cell['zIndex_0']['p_0']['coordinates'][0,:,:]
            
            nrows,ncols = cell_micron_coords.shape
            cell_pixel_coords = np.matmul(
                A_micron_to_mosaic,
                np.vstack((cell_micron_coords.T,np.ones((1,nrows)))),
            ).T
            
            cell_poly = Polygon(cell_pixel_coords[:,:2])
            
            for raw_p in raw_nuclei_polys:
                nuc_pixel_coords = approximate_polygon(raw_p, tolerance=2.5) #simplifying the polygons saves a lot of space and looks good still
                nuc_pixel_coords[:,0] += fov_pixel_min_x
                nuc_pixel_coords[:,1] += fov_pixel_min_y

                #skip polygons without at least 3 vertices
                num_vertices,_ = nuc_pixel_coords.shape
                if num_vertices < 3:
                    continue

                nuc_poly = Polygon(nuc_pixel_coords)
                
                #plotting the good ones
                if cell_poly.contains(nuc_poly):
                    nrows,ncols = nuc_pixel_coords.shape
                    nuc_micron_coords = np.matmul(
                        A_mosaic_to_micron,
                        np.vstack((nuc_pixel_coords.T,np.ones((1,nrows)))),
                    ).T

                    nuc_micron_poly = Polygon(nuc_micron_coords[:,:2])
                    
                    micron_nuclei_polys.append(nuc_micron_poly)
                    cell_ids.append(cell_id)
                    
    gdf = geopandas.GeoDataFrame({
        'cell_id':cell_ids,
        'geometry':micron_nuclei_polys,
    })
    gdf['area'] = gdf.area

    return gdf


def main():
    global im

    #Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample')
    parser.add_argument('--z')
    parser.add_argument('--mosaic')
    parser.add_argument('--transform')
    parser.add_argument('--hdf5_fovs', nargs='+')
    args = parser.parse_args()

    #Read in the DAPI stain image (slowest step?)
    im = io.imread(args.mosaic)

    #Iterate through all fov's of this sample at this z-slice
    #write out in append mode to avoid storing all in memory
    written = False
    for hdf5_fov in args.hdf5_fovs:
        sys.stdout.write(hdf5_fov+'\n')
        sys.stdout.flush()

        fov_gdf = nuclei_segment_fov(args.transform, hdf5_fov)

        if fov_gdf.empty:
            continue

        fov_gdf['sample'] = args.sample
        fov_gdf['z'] = args.z
        fov_gdf['fov'] = os.path.basename(hdf5_fov)
        fov_gdf.to_file('cell_nuc.gpkg', driver='GPKG', mode=('w' if not written else 'a'))
        written = True
 
    #Have to handle no nucleus case specially, not allowed to write an empty file
    #will have to deal with this in the merging script as well
    #(don't actually expect this to happen)
    if not written:
        with open('cell_nuc.gpkg','w') as f:
            f.write('empty file\n')


if __name__ == '__main__':
    main()

