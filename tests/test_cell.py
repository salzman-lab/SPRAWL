import pytest
import SRRS
from SRRS import cell

import shapely.geometry

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
def test_cell_gene_count(dataset, request):
    dataset = request.getfixturevalue(dataset)

    for cell in dataset.iter_cells():
        gc_total = sum(cell.gene_counts.values())
        assert gc_total == cell.n


@pytest.mark.parametrize(
    'dataset,threshold,filt_count', [
        ('m1s4',800,6),
        ('m2s4',800,11),
    ]
)
def test_iter_cells_filter_tot_spots(dataset, threshold, filt_count, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()
    cells = (c for c in cells if c.n > threshold)

    num_filt_cells = sum(1 for c in cells)
    assert num_filt_cells == filt_count


@pytest.mark.parametrize(
    'dataset,threshold,filt_count', [
        ('m1s4',115,9),
        ('m2s4',115,8),
    ]
)
def test_iter_cells_filter_min_unique_genes(dataset, threshold, filt_count, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()
    cells = (c for c in cells if len(c.gene_counts) > threshold)

    num_filt_cells = sum(1 for c in cells)
    assert num_filt_cells == filt_count


@pytest.mark.parametrize(
    'dataset,min_genes,min_spots_per_gene,filt_count', [
        ('m1s4',20,10,14),
        ('m2s4',20,10,18),
    ]
)
def test_iter_cells_filter_min_unique_genes_ge_threshold(dataset, min_genes, min_spots_per_gene, filt_count, request):
    dataset = request.getfixturevalue(dataset)

    cells = dataset.iter_cells()

    #Filtering cells to have:
    #"at least 20 unique genes with each gene having at least 10 RNA spot counts"
    cells = (c for c in cells if sum(v >= min_spots_per_gene for v in c.gene_counts.values()) >= min_genes)

    num_filt_cells = sum(1 for c in cells)
    assert num_filt_cells == filt_count



@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_shrunk_cells_area_decreases(dataset, request):
    dataset = request.getfixturevalue(dataset)

    for cell in dataset.iter_cells():
        shrunk_cell = cell.shrink_boundaries(scale_factor=0.8)

        for zslice in cell.zslices:
            orig_p = shapely.geometry.Polygon(cell.boundaries[zslice])
            shrunk_p = shapely.geometry.Polygon(shrunk_cell.boundaries[zslice])

            assert orig_p.area > shrunk_p.area


