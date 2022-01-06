import SRRS
from SRRS import scoring
import sys

def main():
    hdf5_path = sys.argv[1]

    sample = SRRS.HDF5(hdf5_path)
    cells = sample.iter_cells()

    #Calc gene vars
    cells = scoring._iter_vars(cells)

    #save gene vars
    sample.save_gene_vars(cells)


if __name__ == '__main__':
    main()

