from . import metrics
from . import utils

from statsmodels.stats import multitest
import scipy as scp
import pandas as pd
import numpy as np

import multiprocessing as mp
import functools
import logging
import sys

pd.options.mode.chained_assignment = None  # default='warn'

available_metrics = {
    '_test':metrics._test,
    'peripheral':metrics.peripheral,
    'radial':metrics.radial,
    'punctate':metrics.punctate,
    'central':metrics.central,
}

def iter_scores(cells, metric):
    """
    Apply the chosen scoring metric to each cell
    utilizes multiprocessing

    yields an iterator of dataframes, one for each cell, each with score info
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
    logging.debug('Scoring with {} cpu cores'.format(mp.cpu_count()))

    cells = _iter_vars(cells)
    cells = _iter_scores(cells,metric)

    for cell in cells:
        genes = cell.gene_med_ranks.keys()
        med_ranks = cell.gene_med_ranks.values()

        df = pd.DataFrame({
            'metric':metric,
            'cell_id':cell.cell_id,
            'annotation':cell.annotation,
            'num_spots':cell.n,
            'gene':genes,
            'num_gene_spots':[cell.gene_counts[g] for g in genes],
            'median_rank':med_ranks,
            'score':[utils.score(r,cell.n) for r in med_ranks],
            'variance':[cell.gene_vars[g] for g in genes],
        })
        yield df


def gene_celltype_scoring(srrs_df, min_cells_per_gene_ont=1, extra_cols={}):
    """
    Calculates gene/cell-type SRRS scores

    input in the form of a pd.DataFrame with the following columns (not necessarily in this order)
         'metric' peripheral/polar etc
         'cell_id' the cell id. will be the same for every row
         'annotation' The cell-type annotation
         'num_spots' The number of spots in the whole cell
         'gene' The gene of the row
         'num_gene_spots' The number of spots of this spots gene in the cell
         'median_rank' the observed median rank
         'score' value ranging from -1 to 1
         'variance' theoretical variance under the null

    Output in the form of an aggregated pd.DataFrame where each row is a gene/ont

    Extra columns can be added to the output using extra_cols
    """

    #Table can be passed in as a df or a path to the df
    if isinstance(srrs_df,str):
        srrs_df = pd.read_csv(srrs_df)

    gb_cols = ['gene','annotation']

    #Filter out gene/onts that have too few cells
    srrs_df = srrs_df.groupby(gb_cols).filter(lambda g: g['cell_id'].nunique() >= min_cells_per_gene_ont)

    #calculate the z-scores using Lyapunov approximation
    z_scores = {}
    for k,g in srrs_df.groupby(gb_cols):
        z_score = g['score'].sum()/np.sqrt(g['variance'].sum())
        z_scores[k] = z_score

    srrs_df = srrs_df.set_index(gb_cols)
    srrs_df['z'] = pd.Series(z_scores)
    srrs_df = srrs_df.reset_index()

    #calculate two-sided p-value
    srrs_df['one_sided_p'] = scp.stats.norm.cdf(srrs_df['z'])
    srrs_df['two_sided_p'] = 2*np.minimum(srrs_df['one_sided_p'], 1-srrs_df['one_sided_p'])

    #multiple testing correction
    dedup_df = srrs_df.drop_duplicates(gb_cols)

    #Running into divide by zero errors which I think is the result of an empty table
    if len(dedup_df) > 1:
        passes,adj_p,_,_ = multitest.multipletests(
            dedup_df['two_sided_p'],
            alpha = 0.05, #just need a stand in value
            method = 'fdr_bh',
        )
        dedup_df['bh_corrected_two_sided_p'] = adj_p
        dedup_df = dedup_df.set_index(gb_cols)
        srrs_df = srrs_df.set_index(gb_cols)
        srrs_df['bh_corrected_two_sided_p'] = dedup_df['bh_corrected_two_sided_p']
        srrs_df = srrs_df.reset_index()

    else:
        srrs_df['bh_corrected_two_sided_p'] = None


    #groupby sample/gene/annotation with bh_p
    agg_df = srrs_df.groupby(gb_cols).agg(
        metric = ('metric','first'),
        gene = ('gene','first'),
        annotation = ('annotation','first'),
        num_cells = ('cell_id','nunique'),
        med_gene_spots = ('num_gene_spots','median'),
        med_spots = ('num_spots','median'),
        med_score = ('score','median'),
        z = ('z','first'),
        p = ('two_sided_p','first'),
        bh_p = ('bh_corrected_two_sided_p','first'),
    ).reset_index(drop=True)

    for col,val in extra_cols.items():
        agg_df[col] = val

    return agg_df


def _sequential_iter_scores(cells, metric):
    """
    Apply the chosen scoring metric to each cell
    does NOT use multiprocessing

    yields an iterator of dataframes with score info
    """
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
            'gene':cell.genes,
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
    #If the gene_vars are already calculated, just return the cell
    if cell.gene_vars:
        return cell

    n = cell.n
    for g,m in cell.gene_counts.items():
        if (m,n) not in var_mem:
            var_mem[(m,n)] = utils.calc_var(m,n)
        cell.gene_vars[g] = var_mem[(m,n)]

    return cell


def _iter_vars(cells, processes=2):
    """
    Helper function
    Calculate the theoretical variance of each gene in each cell
    utilizes multiprocessing
    """
    manager = mp.Manager()
    var_mem = manager.dict()

    with mp.Pool(processes=processes) as p:
        f = functools.partial(_cell_var, var_mem = var_mem)
        for result in p.imap_unordered(f, cells):
            yield result


def _iter_scores(cells, metric_name, processes=2):
    """
    Helper function
    Apply the chosen scoring metric to each cell
    utilizes multiprocessing
    """

    #get the metric function or raise KeyError
    metric_f = available_metrics[metric_name]

    #multiplex the scoring
    with mp.Pool(processes=processes) as p:
        for result in p.imap_unordered(metric_f, cells):
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
    NOTE Unused??
    multiprocessing approach to calculating utils.p_two_sided_med

    input is a dict of set like d[(m,n)] = set(med1,med2)

    output is a dict d[(m,n,med)] = p
    """
    with mp.Pool() as p:
        p_twosided_ds = p.map(_calc_p_twosided_helper, m_n_meds.items())
        d = {k: v for d in p_twosided_ds for k, v in d.items()}
        return d

