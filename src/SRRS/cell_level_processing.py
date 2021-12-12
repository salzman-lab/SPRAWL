import multiprocessing as mp
from . import metrics

import logging
import sys

def iter_scores(cells, metric):
    """
    Apply the chosen scoring metric to each cell
    utilizes multiprocessing

    yields an iterator of dataframes with score info
    """
    available_metrics = {
        '_test':metrics._test,
        'peripheral':metrics.peripheral,
    }

    if metric not in available_metrics:
        sys.stderr.write('Metric {} not found in possible metrics {}\n'.format(metric, available_metrics.keys()))
        sys.exit(1)

    logging.info('Scoring with {} cpu cores'.format(mp.cpu_count()))

    with mp.Pool() as p:
        results = p.imap_unordered(available_metrics[metric], cells)
        for result in results:
            yield result


def _sequential_iter_scores(cells, metric):
    """
    Apply the chosen scoring metric to each cell
    does NOT use multiprocessing

    yields an iterator of dataframes with score info
    """
    available_metrics = {
        '_test':metrics._test,
        'peripheral':metrics.peripheral,
    }

    if metric not in available_metrics:
        sys.stderr.write('Metric {} not found in possible metrics {}\n'.format(metric, available_metrics.keys()))
        sys.exit(1)

    for cell in cells:
        result = available_metrics[metric](cell)
        yield result

