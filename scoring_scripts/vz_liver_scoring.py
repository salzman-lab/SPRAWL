import SRRS

import pandas as pd
import time
import sys
import os

def main():
    prog_start = time.time()
    hdf5_path = sys.argv[1]
    out_dir = 'outputs/raw'
    out_name = os.path.basename(hdf5_path).split('.')[0]+'_peripheral.csv'
    out_path = os.path.join(out_dir,out_name)

    sample = SRRS.HDF5(hdf5_path)
    cells = sample.iter_cells()

    #score and write out
    scores = SRRS.iter_scores(cells, metric='peripheral')

    for i,score_df in enumerate(scores):
        score_df.to_csv(out_path, mode=('w' if i==0 else 'a'), index=False, header = i==0)

if __name__ == '__main__':
    main()

