import pytest

from SRRS import vignette

@pytest.fixture
def m1s4():
    return vignette.m1s4_hdf5()

@pytest.fixture
def m2s4():
    return vignette.m2s4_hdf5()

