import pytest

def test_vignette_1(m1s4):
    assert m1s4.num_cells == 20

def test_vignette_2(m2s4):
    assert m2s4.num_cells == 20


def test_vignette_1_diff_2(m1s4, m2s4):
    assert m1s4.cell_ids != m2s4.cell_ids
