import pandas as pd
import sys

import SRRS

samples = [
    ('before_gene_rename/CZB_kidney_111921.hdf5', 'CZB_kidney_111921.hdf5'),
    ('before_gene_rename/CZB_liver.hdf5', 'CZB_liver.hdf5'),
]

ind = int(sys.argv[1])-1
in_path,out_path = samples[ind]

df = pd.read_csv('../barcode_id_to_gene.csv')
df['barcode_id'] = df['barcode_id'].astype(float)
mapping = dict(df[['barcode_id','gene_name']].values)

sys.stdout.write('Starting {}\n'.format(in_path))
sys.stdout.flush()

sample = SRRS.HDF5(in_path)

def rename_genes(c):
    #Replace gene names in all structures
    c.genes = [mapping[float(g)] for g in c.genes]
    
    c.gene_vars = {mapping[float(g)]:v for g,v in c.gene_vars.items()}
    c.gene_med_ranks = {mapping[float(g)]:v for g,v in c.gene_med_ranks.items()}
    c.gene_counts = {mapping[float(g)]:v for g,v in c.gene_counts.items()}
    
    c.spot_genes = {z:[mapping[float(g)] for g in genes] for z,genes in c.spot_genes.items()}
    return c

cells = (rename_genes(c) for c in sample.iter_cells())
SRRS.HDF5.write_cells(cells, out_path)
sys.stdout.write('Finished {}\n'.format(out_path))
sys.stdout.flush()
 


