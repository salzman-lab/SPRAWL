#!/usr/bin/env python3
from rtree import index
import pandas as pd
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Index spots into an RTree')

    parser.add_argument('--spot_path', dest='spot_path', required=True,
                        help='Path to the RNA spot file to be processed')

    args = parser.parse_args()
    return args


# Main #
if __name__ == '__main__':
    args = parse_args()

    #Make a 3D index
    p = index.Property()
    p.dimension = 3
    idx = index.Rtree('rtree',properties=p)

    spots_df_iter = pd.read_csv(args.spot_path, index_col=0, chunksize=10000)

    i = 0
    for spots_df in spots_df_iter:
        for _,p in spots_df.iterrows():
            x = p['global_x']
            y = p['global_y']
            z = p['global_z']
            gene = str(p['gene'])
            idx.insert(i, (x, y, z, x, y, z), obj=gene)
            i += 1

        


