from . import metrics
from . import utils

from statsmodels.stats import multitest
from scipy import stats as scpstats
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

def iter_scores(cells, metric, **kwargs):
    """
    Apply the chosen scoring metric to each cell
    utilizes multiprocessing

    returns a dataframe with score info for
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
    metric_f = available_metrics[metric]
    score_df = metric_f(cells, **kwargs)
    return score_df


def gene_celltype_scoring(srrs_df, gb_cols=['gene','annotation'], min_cells_per_gene_ont=1, extra_cols={}):
    """
    Calculates aggregated SRRS scores at the gene/cell-type level by default
    Can also be calculated at the gene level using gb_cols=['gene']

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
    srrs_df['one_sided_p'] = scpstats.norm.cdf(srrs_df['z'])
    srrs_df['two_sided_p'] = 2*np.minimum(srrs_df['one_sided_p'], 1-srrs_df['one_sided_p'])

    #multiple testing correction
    dedup_df = srrs_df.drop_duplicates(gb_cols)

    #Running into divide by zero errors which I think is the result of an empty table
    if len(dedup_df) > 1:
        passes,adj_p,_,_ = multitest.multipletests(
            dedup_df['two_sided_p'],
            alpha = 0.05,
            method = 'fdr_bh',
        )
        dedup_df['bh_corrected_two_sided_p'] = adj_p
        dedup_df = dedup_df.set_index(gb_cols)
        srrs_df = srrs_df.set_index(gb_cols)
        srrs_df['bh_corrected_two_sided_p'] = dedup_df['bh_corrected_two_sided_p']
        srrs_df = srrs_df.reset_index()

    else:
        srrs_df['bh_corrected_two_sided_p'] = None

    #Add experiment and sample columns if not present
    if 'experiment' not in srrs_df:
        srrs_df['experiment'] = 'experiment1'
    if 'sample' not in srrs_df:
        srrs_df['sample'] = 'sample1'

    #group and dedup with bh_p values
    agg_df = srrs_df.groupby(gb_cols).agg(
        experiment = ('experiment','first'),
        sample = ('sample','first'),
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


