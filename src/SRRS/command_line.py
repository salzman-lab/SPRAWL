import click
from . import scoring

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='0.0.1')
def srrs():
    """
    cli for SRRS to detect subcellular RNA localization patterns
    """
    pass


@srrs.command()
@click.argument('filename', type=click.Path(exists=True))
def validate(filename):
    """
    Validate that [filename] hdf5 is in a usable format for SRRS
    """
    #TODO just echoing back the file name, need to implement the actual logic
    click.echo(click.format_filename(filename))


@srrs.command()
@click.argument('input', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
@click.option('--metric', type=click.Choice(scoring.available_metrics.keys()), required = True)
def gene_cell_scoring(**kwargs):
    """
    Calculate gene/cell scores of an HDF5 [input] file using one of the spatial metrics and save as a CSV to the [output] file
    """
    #TODO not sure if I want to bother adding this functionality to the cli
    input_str = click.format_filename(kwargs['input'])
    output_str = click.format_filename(kwargs['output'])
    click.echo('Will score {} using {} and output to {}'.format(input_str, kwargs['metric'], output_str))



if __name__ == '__main__':
    srrs()

