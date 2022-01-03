from . import metrics
from . import utils

import multiprocessing as mp
import pandas as pd
import numpy as np
import functools
import logging
import sys


def iter_scores(cells, metric):
    """
    Apply the chosen scoring metric to each cell
    utilizes multiprocessing

    yields an iterator of dataframes with score info
     'metric' peripheral/polar etc
     'cell_id' the cell id. will be the same for every row
     'annotation' The cell-type annotation information. same for all rows
     'num_spots' The number of spots in the whole cell. same for all rows
     'gene' The gene of the row
     'num_gene_spots' The number of spots of this spots gene in the cell
     'median_rank' the observed median rank
     'score' value ranging from -1 to 1
     'variance' theoretical variance under the null
    """
    logging.info('Scoring with {} cpu cores'.format(mp.cpu_count()))

    cells = _iter_vars(cells)
    cells = _iter_scores(cells,metric)

    for cell in cells:
        df = pd.DataFrame({
            'metric':metric,
            'cell_id':cell.cell_id,
            'annotation':cell.annotation,
            'num_spots':cell.n,
            'gene':[g.decode() for g in cell.genes],
            'num_gene_spots':[cell.gene_counts[g] for g in cell.genes],
            'median_rank':[cell.gene_med_ranks[g] for g in cell.genes],
            'score':[utils.score(cell.gene_med_ranks[g],cell.n) for g in cell.genes],
            'variance':[cell.gene_vars[g] for g in cell.genes],
        })
        yield df


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
        cell = _cell_var(cell)
        cell = available_metrics[metric](cell)

        df = pd.DataFrame({
            'metric':metric,
            'cell_id':cell.cell_id,
            'annotation':cell.annotation,
            'num_spots':cell.n,
            'gene':[g.decode() for g in cell.genes],
            'num_gene_spots':[cell.gene_counts[g] for g in cell.genes],
            'median_rank':[cell.gene_med_ranks[g] for g in cell.genes],
            'score':[utils.score(cell.gene_med_ranks[g],cell.n) for g in cell.genes],
            'variance':[cell.gene_vars[g] for g in cell.genes],
        })
        yield df


def _cell_var(cell, var_mem={}):
    """
    Helper function to calculate theoretical gene variance
    utilizes manager.dict() shared memory
    adds gene_vars as a member of the cell object
    """
    n = cell.n
    for g,m in cell.gene_counts.items():
        if (m,n) not in var_mem:
            var_mem[(m,n)] = utils.calc_var(m,n)
        cell.gene_vars[g] = var_mem[(m,n)]

    return cell


def _iter_vars(cells):
    """
    Helper function
    Calculate the theoretical variance of each gene in each cell
    utilizes multiprocessing
    """
    manager = mp.Manager()
    var_mem = manager.dict()

    with mp.Pool() as p:
        f = functools.partial(_cell_var, var_mem = var_mem)
        results = p.imap_unordered(f, cells)
        for result in results:
            yield result


def _iter_scores(cells, metric):
    """
    Helper function
    Apply the chosen scoring metric to each cell
    utilizes multiprocessing
    """
    available_metrics = {
        '_test':metrics._test,
        'peripheral':metrics.peripheral,
    }

    #get the metric function or raise KeyError
    metric_f = available_metrics[metric]

    #multiplex the scoring
    with mp.Pool() as p:
        results = p.imap_unordered(metric_f, cells)
        for result in results:
            yield result

def _calc_p_twosided_helper(m_n_meds):
    """
    Helper function to calculate theoretical p_meds

    Input m_n_meds is a tuple of ((m,n),set(med1,med2,med3,...))
    Output is a dict d[(m,n,med)] = two_sided_p
    """
    (m,n),obs_meds = m_n_meds
    obs_meds = list(obs_meds) #avoid order changes in set
    ps = {(m,n,obs_med):0 for obs_med in obs_meds} #init the return dict
    exp_med = (n+1)/2

    #symmetrical, flip to be on < side
    flipped_meds = [n+1-obs_med if obs_med > exp_med else obs_med for obs_med in obs_meds]
    max_med = max(flipped_meds)+1

    if m%2 == 1:
        #m odd case
        min_med = m//2+1
        possible_meds = np.arange(min_med,max_med)

    else:
        #m even case
        min_med = (m+1)/2
        possible_meds = np.arange(min_med,max_med,0.5)

    p = 0
    for med in possible_meds:
        med_p = utils.p_med(m,n,med)
        p += med_p

        if med in flipped_meds:
            ps[(m,n,med)] = 2*p #multiply by two since two-sided
            ps[(m,n,n+1-med)] = 2*p #multiply by two since two-sided

    return ps


def _calc_p_twosided(m_n_meds):
    """
    multiprocessing approach to calculating utils.p_two_sided_med

    input is a dict of set like d[(m,n)] = set(med1,med2)

    output is a dict d[(m,n,med)] = p
    """
    with mp.Pool() as p:
        p_twosided_ds = p.map(_calc_p_twosided_helper, m_n_meds.items())
        d = {k: v for d in p_twosided_ds for k, v in d.items()}
        return d
