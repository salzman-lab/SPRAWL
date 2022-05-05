import pytest
from SRRS import plotting,scoring

import matplotlib

@pytest.mark.parametrize('gene_colors', [{},{'Acta2':'red'}])
@pytest.mark.parametrize('color_by_rank', [True,False])
@pytest.mark.parametrize('metric', [None,'peripheral'])
def test_3D_plot(m1s4,gene_colors,color_by_rank,metric):
    cell = m1s4.cells()[0]
    fig,ax = plotting.plot_cell_3D(cell,gene_colors=gene_colors,color_by_rank=color_by_rank)

    assert type(fig) == matplotlib.figure.Figure


@pytest.mark.parametrize('gene_colors', [{},{'Acta2':'red'}])
@pytest.mark.parametrize('color_by_rank', [True,False])
@pytest.mark.parametrize('metric', [None,'peripheral'])
def test_zslice_plot(m2s4,gene_colors,color_by_rank,metric):
    cell = m2s4.cells()[0]
    fig,axs = plotting.plot_cell_zslices(cell,gene_colors=gene_colors,color_by_rank=color_by_rank)

    assert type(fig) == matplotlib.figure.Figure
    assert len(axs) >= len(cell.zslices) #make sure there is a subplot for each z-slice (sometimes more subplots due to grid)


@pytest.mark.parametrize('color_by_score_gene, color_by_ontology',
        [
            (None,False),
            (None,True),
            ('Acta2',False),
        ]
)
def test_tissue_plot(m2s4,color_by_score_gene,color_by_ontology):
    cells = m2s4.cells()
    fig,ax = plotting.plot_tissue_level(
        cells,
        color_by_score_gene=color_by_score_gene,
        color_by_ontology=color_by_ontology,
    )
    assert type(fig) == matplotlib.figure.Figure


