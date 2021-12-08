import pandas as pd
import numpy as np
import h5py

class HDF5(object):
    """
    Wrapper class around spatial HDF5 file
    """
    def __init__(self, hdf5_path):
        self.hdf5_path = hdf5_path
        self._f = h5py.File(self.hdf5_path,'r')

    def __del__(self):
        if self._f is not None:
            self._f.close()

    @property
    def cell_ids(self):
        return list(self._f['cells'].keys())

    @property
    def num_cells(self):
        return len(self._f['cells'].keys())

    @property
    def genes(self):
        return sorted(g.decode() for g in set(self._f['genes']))

    @property
    def num_genes(self):
        return len(set(self._f['genes']))

