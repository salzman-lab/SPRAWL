import pandas as pd

#read in the gtf here (yao et al use mm10)
full_gtf = pd.read_csv(
    'gencode.vM23.annotation.gtf',
    comment = '#',
    sep = '\t',
    header = None,
    names = ['chrom','source','kind','start','end','dot1','strand','dot2','info'],
)
gtf = full_gtf[full_gtf['kind'].eq('gene')]
gtf['label'] = gtf['info'].str.extract('gene_name "(.*?)";')

#read in the MERFISH genes to allow subsetting
mop_genes = pd.read_csv(
    '/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/preprocessing/BICCN_preprocessing/mop_genes.txt',
    header=None,
    names=['label'],
)
mop_genes = set(mop_genes['label'])

gtf = gtf[gtf['label'].isin(mop_genes)]

#create and writeout the entire-gene bed files (min and max chromosome positions)
entire_gene_bed = gtf.groupby(['label','chrom','strand']).agg(
    start = ('start','min'),
    end = ('end','max'),
).reset_index()
entire_gene_bed['score'] = 0

entire_gene_bed[['chrom','start','end','label','score','strand']].to_csv(
    'BICCN_merf_genes.bed',
    sep=' ',
    header=False,
    index=False,
)

