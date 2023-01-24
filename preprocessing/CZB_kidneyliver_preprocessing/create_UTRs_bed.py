import pandas as pd

#Getting 3' UTR bed file for the Kidney and Liver genes for ReadZs runs


#Read in the mm10 gtf annotated 3' UTR information
gtf_df = pd.read_csv(
    '/oak/stanford/groups/horence/rob/readzs_fork/mm10.refGene.gtf',
    sep='\t',
    header=None,
    names=['chr','source','kind','start','end','dot1','strand','dot2','info'],
)

gtf_df = gtf_df[gtf_df['kind'].eq('3UTR')]
gtf_df['gene'] = gtf_df['info'].str.extract('gene_id "(.*?)";')
gtf_df['gene'] = gtf_df['gene'].str.lower()

gtf_df = gtf_df.groupby(['gene','chr','source','strand']).agg(
    start = ('start','min'),
    end = ('end','max'),
).reset_index()

#Get a list of the gene names used in the CZB MERFISH samples
adata_Kidney = scp.read_h5ad(
    '/scratch/groups/horence/rob/data/Angela_Pisco_MERFISH/cell_gene_counts/MERFISH_kidney_object.h5ad',
)

adata_Liver = scp.read_h5ad(
    '/scratch/groups/horence/rob/data/Angela_Pisco_MERFISH/cell_gene_counts/MERFISH_liver_object.h5ad',
)

merf_genes = set(adata_Kidney.var.index).union(adata_Kidney.var.index)
print(len(merf_genes),'number of genes in CZB MERFISH Liver/Kidney')
merf_genes

#Subset the UTRs
utr_df = gtf_df[gtf_df['gene'].isin(merf_genes)]
print(utr_df.shape[0],'number of found genes in the mm10 gtf')

#create a bed file
utr_df[['chr','start','end','gene','strand']].to_csv('/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/preprocessing/KidneyLiver_preprocessing/CZB_UTRs_mm10.bed',index=False)

