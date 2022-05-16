from statsmodels.stats import multitest
import pandas as pd
import glob

file_stem = '/oak/stanford/groups/horence/rob/readzs_fork/results/medians/20220314_MERF_UTRs_10Xv3_no_unann_filt*'
out_path = '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/outputs/readzs/UTR_level/MOp_10Xv3_no_unann_filt.csv'
bed_path = '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/plotting/buildup_plots/MERFISH_genes.bed'

fps = glob.glob(file_stem)
fps = [fp for fp in fps if 'significant' not in fp]
df = pd.concat(pd.read_csv(fp, sep='\t') for fp in fps)
print('Input df is',df.shape)

#Get mouse and ontology separately
df = df.rename(columns={'ontology':'ontology_mouse'})
df['ontology'] = df['ontology_mouse'].str.split('___').str[0]
df['mouse'] = df['ontology_mouse'].str.split('___').str[1]

#Get genomic window position and gene
split_window = df['window'].str.split('_')
df['chr'] = split_window.str[0]
df['strand'] = split_window.str[2].replace({'plus':'+','minus':'-'})
df['gene'] = split_window.str[1]

#Ignore gene windows from the wrong strand
bed_df = pd.read_csv(bed_path,sep=' ')
expected_strand = df['gene'].map(dict(bed_df[['gene','strand']].values))
print('Before strand filtering',df.shape)
print(df.head())
df = df[df['strand'].eq(expected_strand)]
print('After strand filtering',df.shape)

print('Parsed df is',df.shape)
    
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
    

df.to_csv(out_path,index=False)

