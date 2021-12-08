import pandas as pd
import numpy as np
import scipy as scp
import argparse
import sys

import utils

pd.options.mode.chained_assignment = None  # default='warn'


def quit_early_with_empty_table(out_path):
    """
    Helper function to quit early if there is not data in the df
    Outputs a table with just a header
    """
    columns = [
        'gene',
        'annotation',
        'spot_count_bin',
        'z_score',
        'mean_eff_size',
        'num_cells',
        'min_gene_spots ',
        'med_gene_spots',
        'max_gene_spots',
        'med_tot_spots',
        'med_gene_frac',
        'two_sided_p',
    ]

    pd.DataFrame(columns=columns).to_csv(out_path, index=False)
    sys.exit(0)


################################
#              Main            #
################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform gene-level analysis')

    parser.add_argument('--cell_chunks', dest='cell_chunks', required=True, nargs='+',
                        help='Path(s) to the gene chunks')

    parser.add_argument('--out_path', dest='out_path', required=True,
                        help='Path to CSV for output')

    parser.add_argument('--filtered_cells_out_path', dest='filt_cells_out_path', required=True,
                        help='Path to CSV for filtered cells output')

    parser.add_argument('--min_cells_per_type', dest='cell_type_thresh', default=1, type=int,
                        help='Threshold for cells in a cell-type')

    parser.add_argument('--min_cells_per_type_bin', dest='cell_type_bin_thresh', default=1, type=int,
                        help='Threshold for cells in a cell-type')

    parser.add_argument('--stratify_by_counts', dest='stratify_by_counts', default=True, type=utils.str2bool, const=True, nargs='?',
                        help='Whether or not to stratify by spot count. Default True')

    parser.add_argument('--convert_annotations', dest='convert_anns', default='no', type=str, 
                        help='Convert cell-type annotations. Default "no" but can be "SS2_MOp"')



    args = parser.parse_args()

    print(args)

    #Read in all the cell-chunk files into a single dataframe
    #The ignore_index=True is critical for avoiding overwrite due to: df.loc[g.index,'z_score'] = z_score
    df = pd.concat([pd.read_csv(fpath) for fpath in args.cell_chunks], ignore_index=True)

    #Convert the cell-type annotations if specified
    #TODO only converting using SS2 currently
    if args.convert_anns == 'SS2_MOp':
        df = utils.map_merf_anns_to_seq(df, seq_source='SS2')

    #Filter out cell/gene pairs that are from annots with too few cells
    print('Filtering out cell/gene pairs from annotations with fewer than',args.cell_type_thresh,'cells')
    df = df.groupby(['gene','annotation']).filter(
        lambda g: g['cell_id'].nunique() >= args.cell_type_thresh
    )

    #Check if df is empty, then quit early with a header-only CSV writeout
    if len(df) == 0:
        df.to_csv(args.filt_cells_out_path, index=False)
        quit_early_with_empty_table(args.out_path)

    #Cut into spot count bins
    if args.stratify_by_counts:
        bins = [0,10,20,30,50,100,150] #separate by spot-count
    else:
        bins = [0,df['num_gene_spots'].max()] #lump everything into a single bin

    df['spot_count_bin'] = pd.cut(df['num_gene_spots'], bins)
    df = df.dropna(subset=['spot_count_bin']) #some cells will have more than the top end of counts, drop those
    df['spot_count_bin'] = df['spot_count_bin'].astype(str)

    #Filter out gene/ann/bin groups that have too few cells
    df = df.groupby(['gene','annotation','spot_count_bin']).filter(
        lambda g: g['cell_id'].nunique() >= args.cell_type_bin_thresh
    )

    #Could have an empty df after the previous filter step as well
    if len(df) == 0:
        df.to_csv(args.filt_cells_out_path, index=False)
        quit_early_with_empty_table(args.out_path)

    #Write out filtered cell/gene rows for the cell collection output
    df.to_csv(args.filt_cells_out_path, index=False)
   

    #Iterate through different gene/cell-type/count-bin groups to calculate z-scores
    memo = {}
    for (gene,annotation,count_bin),g in df.groupby(['gene','annotation','spot_count_bin']):

        z_score,memo = utils.normal_mean_effect_sizes(
            df = g,
            m_col='num_gene_spots', 
            n_col='num_spots',
            eff_col='gene_score',
            var_memo=memo,
        )

        df.loc[g.index,'z_score'] = z_score

    #Agg to get a df where each row is a unique gene/cell-type
    df['gene_frac'] = (df['num_gene_spots']/df['num_spots']).astype(float)

    df = df.groupby(['gene','annotation','spot_count_bin','z_score']).agg(
        mean_eff_size = ('gene_score','mean'),
        num_cells = ('cell_id','nunique'),
        min_gene_spots = ('num_gene_spots','min'),
        med_gene_spots = ('num_gene_spots','median'),
        max_gene_spots = ('num_gene_spots','max'),
        med_tot_spots = ('num_spots','median'),
        med_gene_frac = ('gene_frac','median'),
    ).reset_index()

    #Add one-sided and two-sided p-values
    df['one_sided_p'] = df['z_score'].apply(scp.stats.norm.cdf)
    df['two_sided_p'] = 2*np.minimum(df['one_sided_p'], 1-df['one_sided_p'])

    #Write out the result df
    df.to_csv(args.out_path, index=False)



