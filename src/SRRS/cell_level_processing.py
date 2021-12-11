import pandas as pd
import numpy as np

import multiprocessing
from . import metrics

import sys

def score_cells(hdf5, metric, min_gene_spots=10):

    available_metrics = {
        '_test':metrics._test,
        'peripheral':metrics.peripheral,
    }


    if metric not in available_metrics:
        sys.stderr.write('Metric {} not found in possible metrics {}\n'.format(metric, available_metrics.keys()))
        sys.exit(1)

    #h5py objects cannot be passed into multiprocessing since they are not pickleable
    #workaround is for each process to be pased the hdf5 path and the cell name to operate on
    #each process then needs to open the hdf5 file, which is something I'd like to avoid as a future TODO
    with multiprocessing.Pool() as p:
        results = p.map(available_metrics[metric], hdf5.iter_cells())

    return results

