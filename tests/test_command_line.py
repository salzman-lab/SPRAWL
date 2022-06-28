import pytest
import SRRS
from SRRS import vignette, command_line

from click.testing import CliRunner

@pytest.mark.skip(reason='Dont know how to invoke() click subcommand')
@pytest.mark.parametrize('processes',[1,5])
def test_ann_bam(tmp_path,processes):
    runner = CliRunner()

    correct_out_path = vignette.get_data_path('ont_ann.bam')
    bam_path = vignette.get_data_path('no_ont.bam')
    out_path = tmp_path / 'test_ont.bam'
    mapping_path = vignette.get_data_path('CB_XO_mapping.csv')

    result = runner.invoke(command_line.srrs, ['annotate', bam_path, out_path, mapping_path])

    assert result.exit_code == 0


