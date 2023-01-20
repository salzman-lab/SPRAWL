import sys
sys.path.append('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/src') #ugly

import pandas as pd
import SRRS
from SRRS import HDF5

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Rename genes in all cells')

    parser.add_argument('--hdf5_path', dest='hdf5_path', required=True,
                        help='Path to the hdf5 file to be processed')

    parser.add_argument('--out_path', dest='out_path', required=True,
                        help='Path to output hdf5 file')

    parser.add_argument('--mapping_table', dest='mapping_table', required=True,
                        help='CSV formatted table with two columns mapping 1st to 2nd')


    args = parser.parse_args()
    return args



def main():
    args = parse_args()

    #Prepare the mapping dictionary
    d_df = pd.read_csv(args.mapping_table)
    c1,c2 = d_df.columns
    d_df[c1] = d_df[c1].astype(str)
    d_df[c2] = d_df[c2].astype(str)
    d = dict(d_df[[c1,c2]].values)

    #Load the sample, rename all genes in all cells, and write out
    sample = SRRS.HDF5(args.hdf5_path)
    cells = (c.rename_genes(d) for c in sample.iter_cells())
    SRRS.HDF5.write_cells(cells, args.out_path)


main()


