import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from skimage import io
from skimage.segmentation import watershed, mark_boundaries
from skimage.filters import sobel
from skimage.exposure import histogram
from scipy import ndimage as ndi
from skimage.measure import find_contours, approximate_polygon
from shapely.geometry import Point, Polygon
import geopandas


#Following tutorial here
#https://scikit-image.org/docs/stable/user_guide/tutorial_segmentation.html

#read in one of the DAPI files
img = io.imread('/scratch/groups/horence/rob/data/vz_liver_showcase/Liver1Slice1/images/mosaic_DAPI_z0.tif')

#get the transformation matrix
micron_to_mosaic = np.genfromtxt(
    '/scratch/groups/horence/rob/data/vz_liver_showcase/Liver1Slice1/images/micron_to_mosaic_pixel_transform.csv',
    delimiter=' ',
)
mosaic_to_micron = np.linalg.inv(micron_to_mosaic)

#Determine masking thresholds (might have to play with these numbers a lot)
hist, hist_centers = histogram(img)
sns.lineplot(x = hist_centers, y = hist)
plt.title('Histogram of pixel values')
plt.savefig('histogram_of_pixel_values.png')
plt.close()

elevation_map = sobel(img)

markers = np.zeros_like(img)
markers[img < 500] = 1
markers[img > 1500] = 2

segmentation = watershed(elevation_map, markers)
segmentation = ndi.binary_fill_holes(segmentation - 1)

#THIS GIVES BACK (y,x) for some reason!!
#followed tutorial here: https://scikit-image.org/docs/dev/auto_examples/edges/plot_polygon.html#sphx-glr-auto-examples-edges-plot-polygon-py
raw_nuclei_polys = find_contours(segmentation,0)

micron_nuclei_polys = []
for raw_p in raw_nuclei_polys:
    mosaic_coords = approximate_polygon(raw_p, tolerance=2.5) #simplifying the polygons saves a lot of space and looks good still
    
    #skip polygons without at least 3 vertices
    num_vertices,_ = mosaic_coords.shape
    if num_vertices < 3:
        continue
    
    #converting from mosaic to micron coordinates
    nrows,ncols = mosaic_coords.shape
    micron_coords = np.matmul(
        mosaic_to_micron,
        np.vstack((mosaic_coords.T,np.ones((1,nrows)))),
    )
    
    y,x,z = micron_coords
    micron_nuclei_polys.append(Polygon(zip(x,y)))
    

    
#Write out results to a GeoPandas datafram
gdf = geopandas.GeoDataFrame({
    'geometry':micron_nuclei_polys,
})
gdf['area'] = gdf.area
gdf.to_file('L1S1_z0_raw_nuclei.gpkg', driver='GPKG') #GPKG saves out a single file instead of many from shp



