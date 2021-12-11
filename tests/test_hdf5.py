import pytest
import os

import SRRS

def test_cell_counts(m1s4, m2s4):
    assert m1s4.num_cells == 20
    assert m2s4.num_cells == 20


def test_vignettes_differ(m1s4, m2s4):
    assert m1s4.cell_ids != m2s4.cell_ids


def test_annotations(m1s4, m2s4):
    assert 'Pvalb_4' in m1s4.annotations
    assert 'L6_CT_3' in m2s4.annotations


def test_cells(m1s4, m2s4):
    sum(1 for _ in m1s4.iter_cells()) == 20
    sum(1 for _ in m2s4.iter_cells()) == 20

