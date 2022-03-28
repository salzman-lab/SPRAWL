from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
import collections
import numpy as np
import shapely
import shapely.geometry
import copy
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


    def shrink_boundaries(self, scale_factor=0.8):
        """
        Shrink the cell boundaries inwards towards the cell center by a certain factor
        Filter out RNA spots that are now outside the cell
        Return a new cell object
        """
        new_cell = copy.deepcopy(self)

        new_cell.gene_vars = {} #if any spot is lost, gene-vars must be recalculated
        self.spot_ranks = {} #ranking will have to be done new
        self.spot_values = {} #ranking will have to be done new

        self.n_per_z = {}
        new_cell.spot_genes = {}
        new_cell.spot_coords = {}


        for zslice in self.zslices:
            #Shrink the boundary
            s = shapely.geometry.Polygon(self.boundaries[zslice])
            s = shapely.affinity.scale(s, xfact=scale_factor, yfact=scale_factor) #default shrinks around polygon center
            new_cell.boundaries[zslice] = np.array(s.boundary.xy).T

            #Remove cells outside the decreased boundary
            new_cell.n_per_z[zslice] = 0
            xy = []
            gene_list = []

            for gene,(x,y) in zip(self.spot_genes[zslice],self.spot_coords[zslice]):
                p = shapely.geometry.Point(x,y)
                if not s.contains(p):
                    continue

                new_cell.n_per_z[zslice] += 1
                xy.append((x,y))
                gene_list.append(gene)

            new_cell.spot_coords[zslice] = np.array(xy)
            new_cell.spot_genes[zslice] = gene_list
                
            
        #update cell summary info
        new_cell.gene_counts = collections.Counter(g for z in new_cell.zslices for g in new_cell.spot_genes[z])
        new_cell.genes = sorted(list(new_cell.gene_counts.keys()))
        new_cell.gene_med_ranks = {}
        new_cell.n = sum(new_cell.n_per_z.values())


        return new_cell



