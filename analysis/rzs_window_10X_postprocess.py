from statsmodels.stats import multitest
import pandas as pd
import glob

file_stem = '/oak/stanford/groups/horence/rob/readzs_fork/results/medians/20220126_*windows*'
out_path = '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/outputs/readzs/window_level/MOp_10Xv3.csv'
bed_path = '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/plotting/buildup_plots/MERFISH_genes.bed'
window_size = 5000

fps = glob.glob(file_stem)
fps = [fp for fp in fps if 'significant' not in fp]
df = pd.concat(pd.read_csv(fp, sep='\t') for fp in fps)

#Get mouse and ontology separately
df = df.rename(columns={'ontology':'ontology_mouse'})
df['ontology'] = df['ontology_mouse'].str.split('___').str[0]
df['mouse'] = df['ontology_mouse'].str.split('___').str[1]

#Get genomic window position
split_window = df['window'].str.split('_')
df['chr'] = split_window.str[0]
df['strand'] = split_window.str[2].replace({'plus':'+','minus':'-'})
df['bin_num'] = split_window.str[1].astype(int)

#Assign genes to windows
bed_df = pd.read_csv(bed_path,sep=' ')
for i,r in bed_df.iterrows():
    min_bin = r['start']//window_size
    max_bin = r['end']//window_size+1
    gene_inds = (
        df['chr'].eq(r['#chr']) &
        df['strand'].eq(r['strand']) &
        df['bin_num'].between(min_bin,max_bin)
    )
    df.loc[gene_inds,'gene'] = r['gene']
    
    
#BH Multiple hypothesis testing separately for each mouse
p_maps = {}
df['perm_p_val'] = df['perm_p_val'].fillna(1)

for m,g in df.groupby('mouse'):
    alpha = 0.05

    window_df = g.drop_duplicates('window')
    
    
    _,adj_p,_,_ = multitest.multipletests(
        window_df['perm_p_val'],
        alpha = alpha,
        method = 'fdr_bh',
    )
    p_maps[m] = pd.Series(index=window_df['window'],data=adj_p)
    

for m,p_map in p_maps.items():
    df.loc[df['mouse'].eq(m),'rz_bh_corrected'] = df.loc[df['mouse'].eq(m),'window'].map(p_map)
    
#Write the output
df.to_csv(out_path,index=False)

