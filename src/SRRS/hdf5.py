import pandas as pd
import numpy as np
import h5py

from . import cell


class HDF5:
    """
    Wrapper class around spatial HDF5 file

    contains helper functions to access data
    decoupled from HDF5 format
    """
    path = None

    def __init__(self, hdf5_path):
        self.path = hdf5_path

    def __repr__(self):
        str_repr = 'HDF5 {}'.format(self.path)
        return str_repr

    def fhandle(func):
        """
        Decorator to open and close hdf5
        """
        def wrapper(self):
            self._f = h5py.File(self.path,'r')
            ret = func(self)
            self._f.close()
            return ret

        return wrapper

    @property
    @fhandle
    def cell_ids(self):
        """
        List of all cell_ids
        """
        return list(self._f['cells'].keys())

    @property
    @fhandle
    def genes(self):
        """
        List of unique genes in alphabetical order
        """
        return sorted(g.decode() for g in set(self._f['genes']))

    @property
    @fhandle
    def annotations(self):
        """
        List of all annotations of the cells in the HDF5 in order of the cells
        """
        return [cell.attrs['annotation'] for cell_id,cell in self._f['cells'].items()]

    @property
    @fhandle
    def num_cells(self):
        """
        Count the number of cells in the HDF5
        """
        return len(self._f['cells'].keys())

    @property
    @fhandle
    def num_genes(self):
        """
        Count the number of unique genes across all cells
        """
        return len(self.genes)


    def get_cells_by_id(self,*cell_ids):
        """
        Return a single cell object by cell_id
        """
        cells = []
        with h5py.File(self.path,'r') as _f:
            for cell_id in cell_ids:
                if cell_id.encode() not in _f['cell_ids']:
                    return None

                cells.append(cell.Cell(_f['cells'][cell_id]))

        return cells


    def cells(self):
        """
        return a list of Cell objects
        """
        cells = []
        with h5py.File(self.path,'r') as _f:
            for cell_id,cell_group in _f['cells'].items():
                cells.append(cell.Cell(cell_group))

        return cells


    def iter_cells(self):
        """
        produce an iterator of Cell objects
        """
        with h5py.File(self.path,'r') as _f:
            for cell_id,cell_group in _f['cells'].items():
                yield cell.Cell(cell_group)



    def save_gene_vars(self,cells):
        """
        given an iterator of Cell objects
        save the per gene per cell variances to the hdf5 file
        cannot be easily multiplexed

        returns an iterator of cells in case they are needed later
        """
        with h5py.File(self.path,'a') as _f:
            for cell in cells:
                group = _f['cells'][cell.cell_id]

                #delete gene vars dataset if already present
                if 'gene_vars' in group:
                    del group['gene_vars']

                #write out gene vars in gene-name (key) sorted order
                var_list = [v for g,v in sorted(cell.gene_vars.items())]
                group.create_dataset('gene_vars', data=var_list)


        return cells


    @classmethod
    def write_cells(cls,cells,out_path):
        """
        Create an hdf5 file from a list of Cell objects
        saves to out_path

        overwrites if the file already exists
        """
        with h5py.File(out_path,'w') as out_f:
            cell_ids = []
            all_genes = set()

            f_cells = out_f.create_group('cells')
            for cell in cells:
                #update global collections
                cell_ids.append(cell.cell_id)
                all_genes = all_genes.union(cell.genes)

                #create cell and attributes
                f_cell = f_cells.create_group(cell.cell_id)
                f_cell.attrs['annotation'] = cell.annotation
                f_cell.attrs['num_genes'] = len(cell.genes)
                f_cell.attrs['num_spots'] = sum(cell.gene_counts.values())
                f_cell.attrs['zslices'] = cell.zslices

                #create boundaries, sport_coords, and spot_genes
                f_bounds = f_cell.create_group('boundaries')
                f_coords = f_cell.create_group('spot_coords')
                f_genes = f_cell.create_group('spot_genes')

                for z in cell.zslices:
                    f_bounds.create_dataset(z, data = cell.boundaries[z])
                    f_coords.create_dataset(z, data = cell.spot_coords[z])
                    f_genes.create_dataset(z, data = [g.encode() for g in cell.spot_genes[z]])


                #create gene_vars if calculated for the cell (make sure to sort by gene)
                if cell.gene_vars:
                    var_list = [v for g,v in sorted(cell.gene_vars.items())]
                    f_gvars = f_cell.create_dataset('gene_vars', data = var_list)

            out_f.create_dataset('genes', data=[g.encode() for g in sorted(list(all_genes))])
            out_f.create_dataset('cell_ids',data=[cell_id.encode() for cell_id in cell_ids])


