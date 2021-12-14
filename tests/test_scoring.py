import pytest
import SRRS

@pytest.mark.slow
@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['_test','peripheral'])
def test_multiplex_scoring_metrics_m1s4(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)
    score_iter = SRRS.iter_scores(dataset.iter_cells(), metric=metric)
    score_dfs = list(score_iter)

    assert len(score_dfs) == dataset.num_cells

    cell_ids = [c for df in score_dfs for c in df.cell_id.unique()]
    assert sorted(cell_ids) == sorted(dataset.cell_ids)

    genes = list(set(g for df in score_dfs for g in df.gene.unique()))
    assert sorted(genes) == sorted(dataset.genes)

    #for cell in dataset.iter_cells():
    #    cell.
    #score_dfs.set_index(['cell_id','gene'])['num_gene_spots']


@pytest.mark.slow
@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
@pytest.mark.parametrize('metric', ['_test','peripheral'])
def test_sequential_scoring_metrics_m1s4(dataset, metric, request):
    dataset = request.getfixturevalue(dataset)
    score_iter = SRRS._sequential_iter_scores(dataset.iter_cells(), metric=metric)
    score_dfs = list(score_iter)

    assert len(score_dfs) == dataset.num_cells

    cell_ids = [c for df in score_dfs for c in df.cell_id.unique()]
    assert sorted(cell_ids) == sorted(dataset.cell_ids)

    genes = list(set(g for df in score_dfs for g in df.gene.unique()))
    assert sorted(genes) == sorted(dataset.genes)



