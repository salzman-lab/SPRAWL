import SRRS

import pandas as pd
import time
import sys
import os

def main():
    prog_start = time.time()
    hdf5_path = sys.argv[1]
    metrics = ['peripheral','radial','punctate','central']

    for metric in metrics:
        out_stem = os.path.basename(hdf5_path).split('.')[0]+'_{}'.format(metric)

        sample = SRRS.HDF5(hdf5_path)

        #gene/cell score
        scores_iter = SRRS.iter_scores(sample.iter_cells(), metric=metric)
        scores_df = pd.concat(scores_iter, ignore_index=True)
        scores_df.to_csv(out_stem+'_gene_cells.csv', index=False)

        #gene/ont score
        gene_ont_df = gene_celltype_scoring(scores_df)

        #write gene/ont output 
        gene_ont_df.to_csv(out_name, index=False)

if __name__ == '__main__':
    main()

