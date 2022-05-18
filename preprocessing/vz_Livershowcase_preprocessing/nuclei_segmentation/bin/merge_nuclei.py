#!/usr/bin/env python3
import pandas as pd
import geopandas
import argparse
import sys

def main():
    #Parse args
    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_nucs', nargs='+')
    args = parser.parse_args()

    #Merge the geodataframes
    gdf = geopandas.GeoDataFrame()
    for cell_nuc in args.cell_nucs:
        #some files are empty and will cause errors like fiona.errors.DriverError
        #I'll just catch all errors though
        try:
            fov_gdf = geopandas.read_file(cell_nuc)
            gdf = pd.concat((gdf,fov_gdf))
        except:
            continue
                   
    gdf.to_file('nuclei.gpkg', driver='GPKG')


if __name__ == '__main__':
    main()

