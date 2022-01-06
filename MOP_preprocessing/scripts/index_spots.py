#!/usr/bin/env python3
from rtree import index
import pandas as pd
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Index spots into an RTree')

    parser.add_argument('--spot_path', dest='spot_path', required=True,
                        help='Path to the RNA spot file to be processed')

    parser.add_argument('--z', dest='z', required=True, type=int,
                        help='Which z-slice to process')



    args = parser.parse_args()
    return args


# Main #
if __name__ == '__main__':
    args = parse_args()

    z_ind_to_global_z = {
        0:0,
        1:1.5,
        2:3,
        3:4.5,
        4:6,
        5:7.5,
        6:9,
    }

    spots_df = pd.read_csv(args.spot_path, index_col=0)
    global_z = z_ind_to_global_z[args.z]
    spots_df = spots_df[spots_df['global_z'].eq(global_z)]

    idx = index.Rtree('rtree')

    for i,p in spots_df.iterrows():
        x = p['global_x']
        y = p['global_y']
        gene = p['target_molecule_name']
        idx.insert(i, (x, y, x, y), obj=gene)

        


