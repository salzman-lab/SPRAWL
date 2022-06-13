import pytest
from SRRS import plotting,scoring,vignette

import matplotlib
import pandas as pd

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


def test_read_buildup_plot():
    bam_path = vignette.get_data_path('timp3.bam')
    locus = ('chr10',86345473,86349665)

    ann_path = vignette.get_data_path('timp3.gtf')
    ann_df = pd.read_csv(ann_path)

    spatial_path = vignette.get_data_path('timp3_peripheral.csv')
    spatial_df = pd.read_csv(spatial_path)

    fig = plotting.read_buildup_plot(bam_path, locus, ann_df, spatial_df, stratify_tag='XO', min_tag_reads=100)
    assert type(fig) == matplotlib.figure.Figure


