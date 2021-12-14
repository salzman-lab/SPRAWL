import multiprocessing as mp
from . import metrics

import functools
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

    #get the metric function or raise KeyError
    metric_f = available_metrics[metric]

    logging.info('Scoring with {} cpu cores'.format(mp.cpu_count()))

    manager = mp.Manager()
    var_mem = manager.dict()

    with mp.Pool() as p:
        f = functools.partial(metric_f, var_mem = var_mem)
        results = p.imap_unordered(f, cells)
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

