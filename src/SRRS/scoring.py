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

def _calc_p_twosided_helper(m_n_med, p_med_mem={}):
    """
    Helper function to calculate theoretical p_meds
    utilizes manager.dict() shared memory
    returns p_two_sided value for the given m_n_med input tuple
    """
    m,n,obs_med = m_n_med

    exp_med = (n+1)/2

    if (m,n,obs_med) not in p_med_mem:
        p_med_mem[(m,n,obs_med)] = utils.p_med(m,n,obs_med)

    p_eq_med = p_med_mem[(m,n,obs_med)]

    if obs_med > exp_med:
        #symmetrical, flip to be on < side
        obs_med = 2*exp_med-obs_med

    if m%2 == 1:
        #m odd case
        min_med = m//2+1
        possible_meds = np.arange(min_med,obs_med)

    else:
        #m even case
        min_med = (m+1)/2
        possible_meds = np.arange(min_med,obs_med,0.5)

    p_twosided = p_eq_med
    for med in possible_meds:
        if (m,n,med) not in p_med_mem:
            p_med_mem[(m,n,med)] = utils.p_med(m,n,med)

        p_twosided += 2*p_med_mem[(m,n,med)]

    return p_twosided


def _calc_p_twosided(m_n_meds):
    """
    multiprocessing approach to calculating utils.p_two_sided_med

    input is a list of tuples like [(m,n,med), (m,n,med), ... , (m,n,med)]
    """
    manager = mp.Manager()
    p_med_mem = manager.dict()

    with mp.Pool() as p:
        f = functools.partial(_calc_p_twosided_helper, p_med_mem = p_med_mem)
        p_twosideds = p.map(f, m_n_meds) #order is maintained
        return {k:p for k,p in zip(m_n_meds,p_twosideds)}

