import pytest
import SRRS

@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['_test','peripheral'])
def test_scoring_metrics_m1s4(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)
    score_dfs = SRRS.score_cells(dataset, metric=metric)
    assert len(score_dfs) == dataset.num_cells

    cell_ids = [c for df in score_dfs for c in df.cell_id.unique()]
    assert cell_ids == dataset.cell_ids

    genes = list(set(g for df in score_dfs for g in df.gene.unique()))
    assert sorted(genes) == sorted(dataset.genes)



