import SRRS
import pandas as pd
import collections
import multiprocessing as mp

#sample and cell id to cluster
cluster_anns = pd.read_csv(
    '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/analysis/vz_Liver/vz_liver_cluster_assignments_40_neighs_10_pcs.csv',
    index_col = 0,
)
cluster_anns.index.name = 'cell_id'
cluster_anns = cluster_anns.reset_index()
ks = cluster_anns['sample']+cluster_anns['replicate']+cluster_anns['cell_id']
vs = cluster_anns['leiden']

key_to_cluster = {k:v for k,v in zip(ks,vs)}

#barcode id to gene name
bc_gn_df = pd.read_csv(
    '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/preprocessing/vz_Liver_showcase_preprocessing/barcode_id_to_gene_name.txt',
    header = None,
    names = ['barcode_id','gene_name'],
)
bc_gn_df['barcode_id'] = bc_gn_df['barcode_id'].astype(str)

bc_to_gene = dict(bc_gn_df[['barcode_id','gene_name']].values)


samples = [
    ('Liver1','Slice1','finished_outputs/Liver1Slice1.hdf5'),
    ('Liver1','Slice2','finished_outputs/Liver1Slice2.hdf5'),
    ('Liver2','Slice1','finished_outputs/Liver2Slice1.hdf5'),
    ('Liver2','Slice2','finished_outputs/Liver2Slice2.hdf5'),
]

def edit_cells(sample, liv, sli):
    for i,cell in enumerate(sample.iter_cells()):
        #change the annotation
        key = liv+sli+cell.cell_id
        cell.annotation = key_to_cluster[key]
        
        #change the gene names in multiple places
        cell.genes = list(map(bc_to_gene.get, cell.genes))
        cell.spot_genes = {z:list(map(bc_to_gene.get, genes)) for z,genes in cell.spot_genes.items()}
        cell.gene_vars = {bc_to_gene[g]:v for g,v in cell.gene_vars.items()}
        cell.gene_counts = collections.Counter({bc_to_gene[g]:v for g,v in cell.gene_counts.items()})
        
        yield cell

def process_sample(sample_tuple):
    liv,sli,p = sample_tuple
    sample = SRRS.HDF5(p)
    cells = edit_cells(sample, liv, sli)
    SRRS.HDF5.write_cells(cells, p+'.edited')
 
with mp.Pool(4) as p:
    p.map(process_sample, samples)

