import SRRS

import pandas as pd
import sys


df = pd.read_csv('vz_liver_cluster_assignments_40_neighs_10_pcs.csv')
df['sample_replicate'] = df['sample']+df['replicate']

samples = [
    ('Liver1Slice1','../../inputs/hdf5s/vz_Liver_mouse1_slice1.hdf5','vz_Liver_ann_mouse1_slice1.hdf5'),
    ('Liver1Slice2','../../inputs/hdf5s/vz_Liver_mouse1_slice2.hdf5','vz_Liver_ann_mouse1_slice2.hdf5'),
    ('Liver2Slice1','../../inputs/hdf5s/vz_Liver_mouse2_slice1.hdf5','vz_Liver_ann_mouse2_slice1.hdf5'),
    ('Liver2Slice2','../../inputs/hdf5s/vz_Liver_mouse2_slice2.hdf5','vz_Liver_ann_mouse2_slice2.hdf5'),
]

row_num = int(sys.argv[1])-1 #for SLURM_TASK_ARRAY_ID
name,in_path,out_path = samples[row_num]

sys.stdout.write('Starting {}\n'.format(name))
sys.stdout.flush()

sample = SRRS.HDF5(in_path)

ann_map = dict(df[df['sample_replicate'].eq(name)][['Unnamed: 0','leiden']].values)
def annotate(c):
    c.annotation = ann_map[c.cell_id]
    return c

annotated_cells = (annotate(c) for c in sample.iter_cells() if c.cell_id in ann_map)
SRRS.HDF5.write_cells(annotated_cells, out_path)
sys.stdout.write('Finished {}\n'.format(name))
sys.stdout.flush()
    

