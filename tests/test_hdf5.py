import pytest
import SRRS

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

