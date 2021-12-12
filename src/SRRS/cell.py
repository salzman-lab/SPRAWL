from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import collections
import os

class Cell:
    """
    Wrapper around an HDF5 cell group

    copies HDF5 cell data into a simple python object
    reason is to allow for multiprocessing to work on pickleable object
    and allow cell data to "live" outside closed HDF5s
    """

    def __init__(self, group):
        self.annotation = group.attrs.get('annotation')
        self.zslices = group.attrs.get('zslices')
        self.cell_id = os.path.basename(group.name)

        self.tot_spots = 0

        self.boundaries = {}
        self.spot_coords = {}
        self.spot_genes = {}
        for zslice in self.zslices:
            self.boundaries[zslice] = group['boundaries'][zslice][:]
            self.spot_coords[zslice] = group['spot_coords'][zslice][:]
            self.spot_genes[zslice] = group['spot_genes'][zslice][:]


    @property
    def gene_counts(self):
        counts = collections.Counter(g for z in self.zslices for g in self.spot_genes[z])
        return counts

def plot_cell_spots_zslices(hdf5_cell, spot_colors={}, default_spot_color='grey'):

    zslices = hdf5_cell.attrs['zslices']
    fig, axs = plt.subplots(figsize=(8,8),nrows=3,ncols=3,sharex=True,sharey=True)
    axs = axs.flatten()

    global_min_x,global_min_y = None,None
    global_max_x,global_max_y = None,None

    for i,zslice in enumerate(zslices):

        #Add the spots
        colors = []
        for gene in hdf5_cell['spot_genes'][zslice]:
            if gene.decode() in spot_colors:
                colors.append(spot_colors[gene.decode()])
            else:
                colors.append(default_spot_color)

        #Draw the cell outline
        boundary = hdf5_cell['boundaries'][zslice][:]
        min_x,min_y = boundary.min(axis=0)
        max_x,max_y = boundary.max(axis=0)

        if not global_min_x or min_x < global_min_x:
            global_min_x = min_x
        if not global_min_y or min_y < global_min_y:
            global_min_y = min_y
        if not global_max_x or max_x > global_max_x:
            global_max_x = max_x
        if not global_max_y or max_y > global_max_y:
            global_max_y = max_y

        polygon = Polygon(boundary, fill=None)
        axs[i].add_artist(polygon)


        axs[i].scatter(
            x = hdf5_cell['spot_coords'][zslice][:,0],
            y = hdf5_cell['spot_coords'][zslice][:,1],
            alpha = 0.8,
            color = colors,
        )


    for used_i in range(i+1):
        axs[used_i].set_xticks([])
        axs[used_i].set_yticks([])
        axs[used_i].set_xlim(global_min_x,global_max_x)
