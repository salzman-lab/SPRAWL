import pandas as pd
import numpy as np

import multiprocessing
from . import metrics

import sys

def score_cells(cells, metric):

    available_metrics = {
        '_test':metrics._test,
        'peripheral':metrics.peripheral,
    }

    if metric not in available_metrics:
        sys.stderr.write('Metric {} not found in possible metrics {}\n'.format(metric, available_metrics.keys()))
        sys.exit(1)

    with multiprocessing.Pool() as p:
        results = p.map(available_metrics[metric], cells)

    return results

