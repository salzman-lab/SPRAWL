import pytest
import SRRS
from SRRS import scoring

import pandas as pd

@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['_test','peripheral','radial','punctate','central'])
def test_scoring_metrics(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)
    score_df = SRRS.iter_scores(dataset.iter_cells(), metric=metric, num_iterations=10)

    assert sorted(score_df['cell_id'].unique()) == sorted(dataset.cell_ids)


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['peripheral','radial','punctate','central'])
def test_gene_celltype_scoring(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.cells()
    srrs_df = SRRS.iter_scores(cells, metric=metric, num_iterations=10)

    agg_df = scoring.gene_celltype_scoring(srrs_df)

    #Each gene can be present in at most nc cells, the number of input cells
    assert agg_df.groupby('gene')['num_cells'].apply(lambda nc: nc < len(cells)).all()

    #TODO more asserts




