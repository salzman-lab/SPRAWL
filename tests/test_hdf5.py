import pytest

import sprawl
from sprawl import utils

def test_cell_counts(m1s4, m2s4):
    assert m1s4.num_cells == 20
    assert m2s4.num_cells == 20


def test_vignettes_differ(m1s4, m2s4):
    assert m1s4.cell_ids != m2s4.cell_ids


def test_annotations(m1s4, m2s4):
    assert 'Pvalb_4' in m1s4.annotations
    assert 'L6_CT_3' in m2s4.annotations


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_iter_cells(dataset, request):
    dataset = request.getfixturevalue(dataset)

    cell_ids = [c.cell_id for c in dataset.iter_cells()]
    assert cell_ids == dataset.cell_ids

    anns = [c.annotation for c in dataset.iter_cells()]
    assert anns == dataset.annotations

@pytest.mark.slow
@pytest.mark.parametrize('dataset', ['temp_m1s4','temp_m2s4'])
def test_save_gene_vars(dataset, request):
    #check gene vars are not already set
    sample = request.getfixturevalue(dataset)
    cells = sample.cells()
    num_gene_var_entries = sum(len(cell.gene_vars) for cell in cells)
    assert num_gene_var_entries == 0

    #calculate gene vars (slow)
    cells = list(utils._iter_vars(cells))
    pre_cache_vars = {cell.cell_id:cell.gene_vars for cell in cells}
    num_gene_var_entries = sum(len(cell.gene_vars) for cell in cells)
    assert num_gene_var_entries > 0

    #save gene vars
    sample.save_gene_vars(cells)

    #test gene vars are read in
    cells = sample.cells()
    post_cache_vars = {cell.cell_id:cell.gene_vars for cell in cells}
    num_gene_var_entries = sum(len(cell.gene_vars) for cell in cells)
    assert num_gene_var_entries > 0


    #test the vars agree with each other before and after caching
    assert pre_cache_vars.keys() == post_cache_vars.keys()

    for cell_id,pre_d in pre_cache_vars.items():
        post_d = post_cache_vars[cell_id]
        assert pre_d == post_d


@pytest.mark.parametrize('dataset', ['m1s4','m2s4'])
def test_write_cells(dataset, request, tmp_path):
    orig_sample = request.getfixturevalue(dataset)
    orig_cells = orig_sample.cells()
    out_path = tmp_path / 'test_save_out.hdf5'

    sprawl.HDF5.write_cells(orig_cells, out_path)
    saved_sample = sprawl.HDF5(out_path)
    saved_cells = saved_sample.cells()

    assert orig_sample.num_cells == saved_sample.num_cells
    assert [c.cell_id for c in orig_cells] == [c.cell_id for c in saved_cells]
    assert [c.gene_vars for c in orig_cells] == [c.gene_vars for c in saved_cells]

