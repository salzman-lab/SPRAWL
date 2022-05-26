import pytest
import SRRS
from SRRS import vignette, utils

import pandas as pd
import filecmp

@pytest.mark.xfail
def test_ann_bam(tmp_path):

    correct_out_path = vignette.get_data_path('ont_ann.bam')
    bam_path = vignette.get_data_path('no_ont.bam')
    out_path = tmp_path / 'test_ont.bam'

    mapping_path = vignette.get_data_path('CB_XO_mapping.csv')
    mapping_df = pd.read_csv(mapping_path)
    mapping = dict(mapping_df.values)

    utils.map_bam_tag(bam_path, out_path, mapping, processes=1)

    assert filecmp.cmp(out_path, correct_out_path)

