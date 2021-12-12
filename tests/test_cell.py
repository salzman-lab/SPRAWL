import pytest
import SRRS
from SRRS import cell

@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_iter_cells_count(dataset, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()
    num_cells = sum(1 for c in cells)

    assert num_cells == dataset.num_cells


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_iter_cells_type(dataset, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()
    correct_type = all(type(c) == cell.Cell for c in cells)

    assert correct_type


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_iter_cells_filter_tot_spots(dataset, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()
    cells = (c for c in cells if c.tot_spots > 800)

    num_filt_cells = sum(1 for c in cells)
    assert num_filt_cells < dataset.num_cells


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_iter_cells_filter_min_unique_genes(dataset, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()
    cells = (c for c in cells if len(c.gene_counts) > 115)

    num_filt_cells = sum(1 for c in cells)
    assert 0 < num_filt_cells < dataset.num_cells


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_iter_cells_filter_min_unique_genes_ge_threshold(dataset, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()

    #Filtering cells to have:
    #"at least 20 unique genes with each gene having at least 10 RNA spot counts"
    cells = (c for c in cells if sum(v >= 10 for v in c.gene_counts.values()) >= 20)

    num_filt_cells = sum(1 for c in cells)
    assert 0 < num_filt_cells < dataset.num_cells


