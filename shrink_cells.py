import SRRS
from SRRS import scoring

import sys

def main():

    hdf5_path = sys.argv[1]
    out_path = hdf5_path+'.shrunk_0.8'

    sample = SRRS.HDF5(hdf5_path)

    #shrink the cells and save to a new file after also calculating gene-cell variance
    shrunk_cells = (c.shrink_boundaries(scale_factor=0.8) for c in sample.iter_cells())
    shrunk_cells = scoring._iter_vars(shrunk_cells)
    SRRS.HDF5.write_cells(shrunk_cells, out_path)

    #read the HDF5 back in and calculate all metrics
    for metric in ['peripheral','punctate','radial']:
        sample = SRRS.HDF5(out_path)
        out_score_name += '.'+metric

        scores = SRRS.iter_scores(sample.iter_cells(), metric=metric)

        for i,score_df in enumerate(scores):
            score_df.to_csv(out_score_name, mode='a', index=False, header = i==0)




if __name__ == '__main__':
    main()

