import pandas as pd
import sys
import os

import SRRS

ann_map = {
    'Astro_1':'Astro','Astro_2':'Astro','Astro_3':'Astro',
    'Endo':'Endo',
    'L23_IT_1':'L23_IT','L23_IT_2':'L23_IT','L23_IT_3':'L23_IT','L23_IT_4':'L23_IT','L23_IT_5':'L23_IT',
    #'L45_IT_1','L45_IT_2','L45_IT_3','L45_IT_4','L45_IT_5','L45_IT_SSp_1','L45_IT_SSp_2',
    'L56_NP_1':'L56_NP','L56_NP_2':'L56_NP',
    'L5_ET_1':'L5_ET','L5_ET_2':'L5_ET','L5_ET_3':'L5_ET','L5_ET_4':'L5_ET','L5_ET_5':'L5_ET',
    'L5_IT_1':'L5_IT','L5_IT_2':'L5_IT','L5_IT_3':'L5_IT','L5_IT_4':'L5_IT',
    'L6_CT_1':'L6_CT','L6_CT_2':'L6_CT','L6_CT_3':'L6_CT','L6_CT_4':'L6_CT','L6_CT_5':'L6_CT','L6_CT_6':'L6_CT','L6_CT_7':'L6_CT','L6_CT_8':'L6_CT','L6_CT_9':'L6_CT',
    'L6_IT_1':'L6_IT','L6_IT_2':'L6_IT','L6_IT_3':'L6_IT',
    'L6_IT_Car3':'L6_IT_Car3',
    'L6b_1':'L6b','L6b_2':'L6b','L6b_3':'L6b',
    'Lamp5_1':'Lamp5','Lamp5_2':'Lamp5','Lamp5_3':'Lamp5','Lamp5_4':'Lamp5','Lamp5_5':'Lamp5','Lamp5_6':'Lamp5','Lamp5_7':'Lamp5','Lamp5_8':'Lamp5','Lamp5_9':'Lamp5',
    #'Micro_1','Micro_2',
    'OPC':'OPC',
    'Oligo_1':'Oligo','Oligo_2':'Oligo','Oligo_3':'Oligo',
    #'PVM',
    #'Peri',
    'Pvalb_1':'Pvalb','Pvalb_10':'Pvalb','Pvalb_11':'Pvalb','Pvalb_12':'Pvalb','Pvalb_2':'Pvalb','Pvalb_3':'Pvalb','Pvalb_4':'Pvalb','Pvalb_5':'Pvalb','Pvalb_6':'Pvalb','Pvalb_7':'Pvalb','Pvalb_8':'Pvalb','Pvalb_9':'Pvalb',
    'SMC':'SMC',
    'Sncg_1':'Sncg','Sncg_2':'Sncg',
    'Sst_1':'Sst','Sst_2':'Sst','Sst_3':'Sst','Sst_4':'Sst','Sst_5':'Sst','Sst_6':'Sst','Sst_7':'Sst','Sst_8':'Sst','Sst_Chodl':'Sst',
    'VLMC':'VLMC',
    'Vip_1':'Vip','Vip_10':'Vip','Vip_2':'Vip','Vip_3':'Vip','Vip_4':'Vip','Vip_5':'Vip','Vip_6':'Vip','Vip_7':'Vip','Vip_8':'Vip','Vip_9':'Vip',
    #'striatum_1','striatum_2',
    #'unannotated',
    #'ventricle_1','ventricle_2'
}

samples = [
    ('before_reannot/BICCN_mouse1sample1.hdf5', 'BICCN_mouse1sample1.hdf5'),
    ('before_reannot/BICCN_mouse1sample2.hdf5', 'BICCN_mouse1sample2.hdf5'),
    ('before_reannot/BICCN_mouse1sample3.hdf5', 'BICCN_mouse1sample3.hdf5'),
    ('before_reannot/BICCN_mouse1sample4.hdf5', 'BICCN_mouse1sample4.hdf5'),
    ('before_reannot/BICCN_mouse1sample5.hdf5', 'BICCN_mouse1sample5.hdf5'),
    ('before_reannot/BICCN_mouse1sample6.hdf5', 'BICCN_mouse1sample6.hdf5'),
    ('before_reannot/BICCN_mouse2sample1.hdf5', 'BICCN_mouse2sample1.hdf5'),
    ('before_reannot/BICCN_mouse2sample2.hdf5', 'BICCN_mouse2sample2.hdf5'),
    ('before_reannot/BICCN_mouse2sample3.hdf5', 'BICCN_mouse2sample3.hdf5'),
    ('before_reannot/BICCN_mouse2sample4.hdf5', 'BICCN_mouse2sample4.hdf5'),
    ('before_reannot/BICCN_mouse2sample5.hdf5', 'BICCN_mouse2sample5.hdf5'),
    ('before_reannot/BICCN_mouse2sample6.hdf5', 'BICCN_mouse2sample6.hdf5'),
]

ind = int(sys.argv[1])-1
in_path,out_path = samples[ind]


sys.stdout.write('Starting {}\n'.format(in_path))
sys.stdout.flush()

sample = SRRS.HDF5(in_path)

def annotate(c):
    c.annotation = ann_map[c.annotation]
    return c

annotated_cells = (annotate(c) for c in sample.iter_cells() if c.annotation in ann_map)
SRRS.HDF5.write_cells(annotated_cells, out_path)
sys.stdout.write('Finished {}\n'.format(out_path))
sys.stdout.flush()
 


