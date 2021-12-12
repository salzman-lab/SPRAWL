import multiprocessing as mp
from . import metrics

import logging
import sys

def iter_scores(cells, metric):

    available_metrics = {
        '_test':metrics._test,
        'peripheral':metrics.peripheral,
    }

    if metric not in available_metrics:
        sys.stderr.write('Metric {} not found in possible metrics {}\n'.format(metric, available_metrics.keys()))
        sys.exit(1)

    logging.info('Scoring with {} cpu cores'.format(mp.cpu_count()))

    with mp.Pool() as p:
        results = p.imap(available_metrics[metric], cells)
        for result in results:
            yield result

