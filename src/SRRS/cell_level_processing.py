import pandas as pd
import numpy as np
import argparse
import h5py

import metrics
import utils

################################
#              Main            #
################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ranks per spot')

    parser.add_argument('--hdf5_path', dest='hdf5_path', required=True,
                        help='Path to the hdf5 file to be processed')

    parser.add_argument('--metric', dest='metric', required=True,
                        choices=metrics.metric_lookup.keys(),
                        help='Which spatial metric to calculate')

    parser.add_argument('--out_stem', dest='out_stem', required=True,
                        help='Path to CSV for output')

    parser.add_argument('--chunk', dest='chunk', default=1, type=int,
                        help='Current cell chunk number to process')

    parser.add_argument('--num_chunks', dest='num_chunks', default=1, type=int,
                        help='Total cell chunks to process')

    parser.add_argument('--gene_chunks', dest='gene_chunks', default=1, type=int,
                        help='Total gene chunks to output into')

    parser.add_argument('--min_gene_spots', dest='min_gene_spots', default=1, type=int,
                        help='Threshold for minimum spots per gene/cell')


    args = parser.parse_args()

    print(args)

    with h5py.File(args.hdf5_path,'r') as f:
        #Determine which cell-ids to process based on the chunk params
        all_cell_ids = f['cell_ids']
        chunk_size = len(all_cell_ids)//args.num_chunks
        chunk_size += 0 if len(all_cell_ids)%args.num_chunks==0 else 1 #account for leftovers        

        chunk_ind_start = (args.chunk-1)*chunk_size
        chunk_ind_end = args.chunk*chunk_size
        cell_ids = all_cell_ids[chunk_ind_start:chunk_ind_end]

        #Also capture which gene groups to process based on the num_gene_chunks
        all_genes = f['genes']
        nc = args.gene_chunks
        cs = len(all_genes)//nc
        cs += 0 if len(all_genes)%nc==0 else 1 #account for leftovers
        gene_to_chunk = {v.decode():i//cs for i,v in enumerate(all_genes)}
        first_write = {i:True for i in range(args.gene_chunks)}


    #Calculate the per-gene score of each cell for the given metric
    #This is an iterator that returns a df for each cell
    per_cell_gene_scores_iter = metrics.calculate_ranks(
        hdf5_path = args.hdf5_path,
        metric = args.metric,
        cell_ids = cell_ids,
    )

    #Group the cells into chunks to process multiple at a time
    #should be a speedup but haven't explicitly checked
    cells_per_chunk = 50
    chunk_cell_gene_scores_iter = utils.grouper(
        per_cell_gene_scores_iter,
        n=cells_per_chunk,
        fillvalue=pd.DataFrame(),
    )

    #Keep track of whether each file has been written to already
    first_write = {i:True for i in range(args.gene_chunks)}

    for gene_scores_iter in chunk_cell_gene_scores_iter:
        #Concat all the chunks together
        gene_scores_df = pd.concat(gene_scores_iter)

        #Filter out cell/gene pairs that have too few spots
        gene_scores_df = gene_scores_df[
            gene_scores_df['num_gene_spots'].ge(args.min_gene_spots)
        ]

        #Filter out cells with too few total spots now
        #TODO

        #Write gene parts to the correct output file
        gene_scores_df['gene_chunk'] = gene_scores_df['gene'].map(gene_to_chunk)

        for chunk,g in gene_scores_df.groupby('gene_chunk'):
            out_path = args.out_stem+str(chunk)

            g.to_csv(
                out_path,
                mode='a',
                header=first_write[chunk],
                index=False,
            )
            first_write[chunk] = False
        

    #Write headers out to the files that didn't have any results
    #This ensures that the transpose operation in nextflow groups the right files
    for chunk,unwritten in first_write.items():
        if not unwritten:
            continue

        header = [
            'gene','gene_score','mn_var','annotation',
            'subclass','num_genes','num_gene_spots',
            'num_spots','volume','cell_id','metric','gene_chunk',
        ]

        out_path = args.out_stem+str(chunk)
        with open(out_path,'w') as f_out:
            f_out.write(','.join(header)+'\n')
        

