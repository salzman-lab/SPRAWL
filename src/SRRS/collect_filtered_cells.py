import pandas as pd
import argparse
import sys


################################
#              Main            #
################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Post processing script')

    parser.add_argument('--gene_chunks', dest='gene_chunks', required=True, nargs='+',
                        help='Path(s) to the gene chunks')

    parser.add_argument('--out_path', dest='out_path', required=True,
                        help='Path to CSV for output')

    args = parser.parse_args()


    #Read in all the different gene files
    df = pd.concat(pd.read_csv(fpath) for fpath in args.gene_chunks)

    #Write out file
    df.to_csv(args.out_path, index=False)
    

