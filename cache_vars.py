import SRRS
from SRRS import scoring
import sys

def main():
    hdf5_path = 'srs/SRRS/vignette_data/merfish_m2s4_vignette.hdf5'

    sample = SRRS.HDF5(hdf5_path)
    cells = sample.iter_cells()

    #Calc gene vars
    cells = scoring._iter_vars(cells)

    #save gene vars
    sample.save_gene_vars(cells)


if __name__ == '__main__':
    main()


