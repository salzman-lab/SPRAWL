import pytest
import os
import shutil

from sprawl import vignette, hdf5

@pytest.fixture
def m1s4():
    return vignette.m1s4_hdf5()

@pytest.fixture
def m2s4():
    return vignette.m2s4_hdf5()

@pytest.fixture
def no_var_m1s4():
    return vignette.no_var_m1s4_hdf5()

@pytest.fixture
def no_var_m2s4():
    return vignette.no_var_m2s4_hdf5()


@pytest.fixture
def temp_m1s4(no_var_m1s4,tmp_path):
    name = os.path.basename(no_var_m1s4.path)
    shutil.copyfile(no_var_m1s4.path, tmp_path / name)
    return hdf5.HDF5(tmp_path / name)

@pytest.fixture
def temp_m2s4(no_var_m2s4,tmp_path):
    name = os.path.basename(no_var_m2s4.path)
    shutil.copyfile(no_var_m2s4.path, tmp_path / name)
    return hdf5.HDF5(tmp_path / name)

