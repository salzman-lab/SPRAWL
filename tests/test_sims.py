import pytest
from SRRS import simulate
import collections
import numpy as np

@pytest.mark.parametrize('dataset',['m1s4','m2s4'])
@pytest.mark.parametrize('seed',[1,77])
def test_perm_by_zslice(dataset,seed,request):
    np.random.seed(seed) #<- want to avoid random rare case where permuting doesn't change gene labels

    data = request.getfixturevalue(dataset)
    orig_cells = data.iter_cells()
    perm_cells = data.iter_cells()

    for orig_c,perm_c in zip(orig_cells,perm_cells):
        #permute the perm_c in place
        simulate.null_permute_gene_labels(perm_c, within_z=True)
        assert list(orig_c.zslices) == list(perm_c.zslices)

        all_same = True
        for z in orig_c.zslices:
            orig_g = orig_c.spot_genes[z]
            perm_g = perm_c.spot_genes[z]

            #same gene label counts per slice before/after permute
            assert collections.Counter(orig_g) == collections.Counter(perm_g)

            #keeping track of whether any permuting occurred
            all_same &= (orig_g == perm_g).all()

        assert not all_same


@pytest.mark.parametrize('dataset',['m1s4','m2s4'])
@pytest.mark.parametrize('seed',[25,682])
def test_perm_across_zslice(dataset,seed,request):
    np.random.seed(seed) #<- want to avoid random rare case where permuting doesn't change gene labels

    data = request.getfixturevalue(dataset)
    orig_cells = data.iter_cells()
    perm_cells = data.iter_cells()

    for orig_c,perm_c in zip(orig_cells,perm_cells):
        #permute the perm_c in place
        simulate.null_permute_gene_labels(perm_c, within_z=True)
        assert list(orig_c.zslices) == list(perm_c.zslices)

        all_same = True
        perm_across_z = False
        for z in orig_c.zslices:
            orig_g = orig_c.spot_genes[z]
            perm_g = perm_c.spot_genes[z]

            #make sure the same number of gene labels are in each slice
            assert len(orig_g) == len(perm_g)

            #keeping track of whether any permuting occurred
            all_same &= (orig_g == perm_g).all()

        assert not all_same


def test_sim_null(m1s4):
    cells = m1s4.iter_cells()
    result = simulate.sim_null_peripheral(cells, within_z=True, n_its=10, alpha=0.05)

    assert set(result['cell_id'].values) == set(m1s4.cell_ids)


