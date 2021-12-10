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


def test_scoring_test(m1s4, m2s4):
    m1s4_score_dfs = SRRS.score_cells(m1s4, metric='_test')
    m2s4_score_dfs = SRRS.score_cells(m2s4, metric='_test')

    assert len(m1s4_score_dfs) == m1s4.num_cells
    assert len(m2s4_score_dfs) == m2s4.num_cells

    m1s4_cell_ids = [c for df in m1s4_score_dfs for c in df.cell_id.unique()]
    m2s4_cell_ids = [c for df in m2s4_score_dfs for c in df.cell_id.unique()]

    assert m1s4_cell_ids == m1s4.cell_ids
    assert m2s4_cell_ids == m2s4.cell_ids

    m1s4_genes = [c for df in m1s4_score_dfs for c in df.gene.unique()]
    m2s4_genes = [c for df in m2s4_score_dfs for c in df.gene.unique()]

    assert sorted(m1s4.genes) == sorted(m1s4.genes)
    assert sorted(m2s4.genes) == sorted(m2s4.genes)


def test_scoring_peripheral(m1s4, m2s4):
    m1s4_score_dfs = SRRS.score_cells(m1s4, metric='peripheral')
    m2s4_score_dfs = SRRS.score_cells(m2s4, metric='peripheral')

    assert len(m1s4_score_dfs) == m1s4.num_cells
    assert len(m2s4_score_dfs) == m2s4.num_cells

    m1s4_cell_ids = [c for df in m1s4_score_dfs for c in df.cell_id.unique()]
    m2s4_cell_ids = [c for df in m2s4_score_dfs for c in df.cell_id.unique()]

    assert m1s4_cell_ids == m1s4.cell_ids
    assert m2s4_cell_ids == m2s4.cell_ids

    m1s4_genes = [c for df in m1s4_score_dfs for c in df.gene.unique()]
    m2s4_genes = [c for df in m2s4_score_dfs for c in df.gene.unique()]

    assert sorted(m1s4.genes) == sorted(m1s4.genes)
    assert sorted(m2s4.genes) == sorted(m2s4.genes)



