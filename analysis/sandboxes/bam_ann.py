import SRRS
from SRRS import utils

import pandas as pd
import sys

f_name = sys.argv[1]

bam_path = '/scratch/groups/horence/rob/data/MERFISH_scRNAseq/10X_mapping/{}.bam'.format(f_name)
out_path = 'ont_{}.bam'.format(f_name)

#Get mapping from cell_bc to celltype
df = pd.read_csv(
    '/oak/stanford/groups/horence/rob/readzs_fork/MOp_10Xv3_metadata.tsv',
    sep = '\t',
)
df = df[df['library'].eq(f_name)]
df['cell_bc'] = df['cell_bc']+'-1'
df['subclass_label'] = df['subclass_label'].map(lambda s: s.replace(' ','_').replace('/',''))

mapping = dict(df[['cell_bc','subclass_label']].values)

utils.map_bam_tag(bam_path, out_path, mapping, processes=12)

