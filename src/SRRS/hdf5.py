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


    def iter_cells(self):
        """
        produce an iterator of Cell objects
        """
        with h5py.File(self.path,'r') as _f:
            for cell_id,cell in _f['cells'].items():
                yield cell.Cell(cell)

