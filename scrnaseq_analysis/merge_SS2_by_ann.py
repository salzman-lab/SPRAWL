#!/bin/bash
#############################
# File Name : merge_SS2_by_ann.py
#
# Purpose : [???]
#
# Creation Date : 24-11-2021
#
# Last Modified : Wed 24 Nov 2021 05:21:15 PM PST
#
# Created By : Rob Bierman
#
##############################

import pandas as pd
import subprocess
import pysam
import glob
import os

#SS2 merging bams by cell-type

#Get list of mapped bam files (bam.bai)
mapped_ss2_paths = glob.glob('/oak/stanford/groups/horence/rob/MERFISH_scRNAseq/ss2_bam_aligns/*.bai')
mapped_ss2_paths = [p.replace('.bai','') for p in mapped_ss2_paths]
ss2_cell_ids = [os.path.basename(p).replace('.bam','') for p in mapped_ss2_paths]

bam_path_df = pd.DataFrame({
    'cell_id':ss2_cell_ids,
    'file_path':mapped_ss2_paths,
})

#Prepare metadata
ss2_ann_df = pd.read_csv(
    '/oak/stanford/groups/horence/rob/MERFISH_scRNAseq/SS2.cluster.annotation.csv',
)
ss2_mem_df = pd.read_csv(
    '/oak/stanford/groups/horence/rob/MERFISH_scRNAseq/SS2.cluster.membership.csv',
    header = None,
    skiprows = 1,
    names = ['cell_id', 'cluster_id'],
)
ss2_meta_df = ss2_mem_df.merge(ss2_ann_df)


#Subset metadata to mapped samples
bam_path_df = bam_path_df.merge(ss2_meta_df, how='left').dropna()


#For each each subclass_label group create a merged bam and create a job
m_cmds = []
i_cmds = []
for ann_name,g in bam_path_df.groupby('subclass_label'):
    print(ann_name)
    out_path = '/oak/stanford/groups/horence/rob/MERFISH_scRNAseq/SS2_merged_bams/{}.bam'.format(ann_name)
    
    #Merge the bams (maintains sorted) and index with a single command
    
    m_cmd = ['samtools','merge','-f','-@','4',out_path]
    m_cmd += g['file_path'].to_list()
    m_cmds.append(m_cmd)
    
    i_cmds.append(['samtools','index',out_path])
    
#Start merging jobs
m_ps = []
for m_cmd in m_cmds:
    p = subprocess.Popen(m_cmd)
    m_ps.append(p)
    
#Wait on the merge jobs to finish
[p.communicate() for p in m_ps]
    
#Start indexing jobs
i_ps = []
for i_cmd in i_cmds:
    p = subprocess.Popen(i_cmd)
    i_ps.append(p)
    
#Wait on the index jobs to finish
[p.communicate() for p in i_ps]

#done!

