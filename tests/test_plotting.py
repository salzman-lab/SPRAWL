import pytest
from SRRS import plotting,scoring

import matplotlib

@pytest.mark.parametrize('gene_colors', [{},{'Acta2':'red'}])
@pytest.mark.parametrize('color_by_rank', [True,False])
@pytest.mark.parametrize('metric', [None,'peripheral'])
def test_3D_plot(m1s4,gene_colors,color_by_rank,metric):
    cell = m1s4.cells()[0]
    if metric:
        #coloring by metric requires ranks to be calculated (defaults to grey otherwise)
        metric_f = scoring.available_metrics[metric]
        metric_f(cell)

    fig = plotting.plot_cell_3D(cell,gene_colors=gene_colors,color_by_rank=color_by_rank)
    assert type(fig) == matplotlib.figure.Figure


@pytest.mark.parametrize('gene_colors', [{},{'Acta2':'red'}])
@pytest.mark.parametrize('color_by_rank', [True,False])
@pytest.mark.parametrize('metric', [None,'peripheral'])
def test_zslice_plot(m2s4,gene_colors,color_by_rank,metric):
    cell = m2s4.cells()[0]
    if metric:
        metric_f = scoring.available_metrics[metric]
        metric_f(cell)

    fig = plotting.plot_cell_zslices(cell,gene_colors=gene_colors,color_by_rank=color_by_rank)
    assert type(fig) == matplotlib.figure.Figure


