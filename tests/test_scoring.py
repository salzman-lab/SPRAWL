import pytest
import SRRS
from SRRS import scoring

import pandas as pd

@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['_test','peripheral',
    pytest.param('radial', marks=pytest.mark.xfail(reason='Radial not implemented'))])
def test_multiplex_scoring_metrics_m1s4(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)
    score_iter = SRRS.iter_scores(dataset.iter_cells(), metric=metric)
    score_dfs = list(score_iter)

    assert len(score_dfs) == dataset.num_cells

    cell_ids = [c for df in score_dfs for c in df.cell_id.unique()]
    assert sorted(cell_ids) == sorted(dataset.cell_ids)

    genes = list(set(g for df in score_dfs for g in df.gene.unique()))
    assert sorted(genes) == sorted(dataset.genes)


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['_test','peripheral',
    pytest.param('radial', marks=pytest.mark.xfail(reason='Radial not implemented'))])
def test_sequential_scoring_metrics_m1s4(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)
    score_iter = SRRS._sequential_iter_scores(dataset.iter_cells(), metric=metric)
    score_dfs = list(score_iter)

    assert len(score_dfs) == dataset.num_cells

    cell_ids = [c for df in score_dfs for c in df.cell_id.unique()]
    assert sorted(cell_ids) == sorted(dataset.cell_ids)

    genes = list(set(g for df in score_dfs for g in df.gene.unique()))
    assert sorted(genes) == sorted(dataset.genes)


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['peripheral',
    pytest.param('radial', marks=pytest.mark.xfail(reason='Radial not implemented'))])
def test_gene_celltype_scoring(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.cells()
    score_iter = SRRS.iter_scores(cells, metric=metric)
    srrs_df = pd.concat(score_iter)

    agg_df = scoring.gene_celltype_scoring(srrs_df)

    #Each gene can be present in at most nc cells, the number of input cells
    assert agg_df.groupby('gene')['num_cells'].apply(lambda nc: nc < len(cells)).all()

    #TODO more asserts




