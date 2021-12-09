import pandas as pd
import numpy as np
import h5py

class HDF5(object):
    """
    Wrapper class around spatial HDF5 file
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
        return list(self._f['cells'].keys())

    @property
    @fhandle
    def genes(self):
        return sorted(g.decode() for g in set(self._f['genes']))

    @property
    @fhandle
    def annotations(self):
        return [cell.attrs['annotation'] for cell_id,cell in self._f['cells'].items()]

    @property
    @fhandle
    def num_cells(self):
        return len(self._f['cells'].keys())

    @property
    @fhandle
    def num_genes(self):
        return len(set(self._f['genes']))



