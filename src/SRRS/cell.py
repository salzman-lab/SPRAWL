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

        self.n = 0
        self.n_per_z = {}
        self.ranked = False
        self.spot_ranks = {}
        self.spot_values = {}

        #boundaries, spot_coords, and spot_gene groups will always be present
        self.boundaries = {}
        self.spot_coords = {}
        self.spot_genes = {}
        for zslice in self.zslices:
            self.n_per_z[zslice] = group['spot_genes'][zslice].shape[0]

            self.boundaries[zslice] = group['boundaries'][zslice][:]
            self.spot_coords[zslice] = group['spot_coords'][zslice][:]
            self.spot_genes[zslice] = [g.decode() for g in group['spot_genes'][zslice][:]]
            self.spot_ranks[zslice] = []
            self.spot_values[zslice] = []

            self.n += self.n_per_z[zslice]

        #members used by other functions
        self.gene_counts = collections.Counter(g for z in self.zslices for g in self.spot_genes[z])
        self.genes = sorted(list(self.gene_counts.keys()))
        self.gene_med_ranks = {}

        #theoretical median gene_vars might be present if the calculation has been done prior
        #gene_vars will correspond with unique genes in the cell in alphabetical order
        self.gene_vars = {}
        if 'gene_vars' in group:
            for i,gene_var in enumerate(group['gene_vars']):
                self.gene_vars[self.genes[i]] = gene_var


    def __repr__(self):
        repr_str = 'Cell-{}-{}'.format(self.cell_id,self.annotation)
        return repr_str

