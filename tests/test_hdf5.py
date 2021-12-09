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

def test_cell_scoring_test(m1s4, m2s4):
    m1s4_scores = SRRS.score_cells(m1s4, metric='_test')
    m2s4_scores = SRRS.score_cells(m2s4, metric='_test')

    assert m1s4_scores == m1s4.cell_ids
    assert m2s4_scores == m2s4.cell_ids


def test_cell_scoring_slow_test(m1s4, m2s4):
    m1s4_scores = SRRS.score_cells(m1s4, metric='_slow_test')
    m2s4_scores = SRRS.score_cells(m2s4, metric='_slow_test')

    assert m1s4_scores == m1s4.cell_ids
    assert m2s4_scores == m2s4.cell_ids


