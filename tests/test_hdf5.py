import pytest

from SRRS import scoring

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


@pytest.mark.parametrize('dataset', ['temp_m1s4','temp_m2s4'])
def test_save_gene_vars(dataset, request):
    #check gene vars are not already set
    dataset = request.getfixturevalue(dataset)
    cells = list(dataset.iter_cells())
    num_gene_var_entries = sum(len(cell.gene_vars) for cell in cells)
    assert num_gene_var_entries == 0
            
    #calculate gene vars (slow)
    cells = list(scoring._iter_vars(cells))
    num_gene_var_entries = sum(len(cell.gene_vars) for cell in cells)
    assert num_gene_var_entries > 0
 
    #save gene vars
    dataset.save_gene_vars(cells)

    #test gene vars are read in
    cells = list(dataset.iter_cells())
    num_gene_var_entries = sum(len(cell.gene_vars) for cell in cells)
    assert num_gene_var_entries > 0 
 

