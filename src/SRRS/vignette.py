import importlib_resources as pkg_resources

from .hdf5 import HDF5
from . import vignette_data

def _load_hdf5(name):
    with pkg_resources.files(vignette_data).joinpath(name) as p:
        hdf5 = HDF5(p)
    return hdf5

def m1s4_hdf5():
    return _load_hdf5('merfish_m1s4_vignette.hdf5')

def m2s4_hdf5():
    return _load_hdf5('merfish_m2s4_vignette.hdf5')

def no_var_m1s4_hdf5():
    #var is not cached
    return _load_hdf5('no_var_merfish_m1s4_vignette.hdf5')

def no_var_m2s4_hdf5():
    #var is not cached
    return _load_hdf5('no_var_merfish_m2s4_vignette.hdf5')

