import pytest
import SRRS
from SRRS import vignette, utils

import pysam
import pandas as pd

def test_ann_bam(tmp_path):
    correct_out_path = vignette.get_data_path('ont_ann.bam')
    bam_path = vignette.get_data_path('no_ont.bam')
    out_path = tmp_path / 'test_ont.bam'

    mapping_path = vignette.get_data_path('CB_XO_mapping.csv')
    mapping_df = pd.read_csv(mapping_path)
    mapping = dict(mapping_df.values)

    utils.map_bam_tag(bam_path, out_path, mapping, processes=1)
    
    with pysam.AlignmentFile(correct_out_path) as true_bam, pysam.AlignmentFile(out_path) as test_bam:
        r1s,r2s = true_bam.fetch(),test_bam.fetch()
        
        while True:
            r1,r2 = next(r1s,False),next(r2s,False)
            
            if not r1 and not r2:
                #both were the same length and hadn't had any mismatches
                assert True
                break

            elif r1 and r2:
                matching = all((
                    r1.seq == r2.seq,
                    r1.query_name == r2.query_name,
                    r1.mapq == r2.mapq,
                    r1.get_tags() == r2.get_tags()
                ))
                if not matching:
                    #some sort of matching issue
                    assert 'Read pair' == 'inconsistency' 
                    break

            else:
                #they were not the same length!
                assert 'Read length' == 'mismatch'
                break

@pytest.mark.parametrize(
    'n,items,result', [
        (1,{'a':10,'b':90,'c':50,'d':40},{('a','b','c','d'):190}),
        (2,{'a':10,'b':90,'c':50,'d':40},{('a','c','d'):100,('b',):90}),
        (3,{'a':10,'b':90,'c':50,'d':40},{('a','d'):50,('c',):50,('b',):90}),
        (4,{'a':10,'b':90,'c':50,'d':40},{('a',):10,('d',):40,('c',):50,('b',):90}),
        (20,{'a':10,'b':90,'c':50,'d':40},{('a',):10,('d',):40,('c',):50,('b',):90}),
    ]
)
def test_balanced_groupings(n,items,result):
    assert utils.create_balanced_groupings(n,items) == result

