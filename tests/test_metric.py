import pytest
import sprawl
from sprawl import scoring, vignette

import pandas as pd

@pytest.mark.xfail
def test_radial():
    #made a fake cell with a known radial pattern for fake gene g_0
    sample = vignette.radial_test_hdf5()
    cells = sample.cells()

    score_df = sprawl.iter_scores(cells, metric='radial')

    assert type(score_df) == pd.core.frame.DataFrame
    assert score_df.sort_values('score',ascending=False).iloc[0]['gene'] == 'g_0' #know that gene0 is the most radial


