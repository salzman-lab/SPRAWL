#!/usr/local/bin/python
#Above is python path in the docker container with SRRS installed
import SRRS
import argparse

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import seaborn as sns

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene_cell_path')
    parser.add_argument('--metric')
    parser.add_argument('--experiment')
    parser.add_argument('--output_name')
    args = parser.parse_args()

    df = pd.read_csv(args.gene_cell_path)

    pdf = matplotlib.backends.backend_pdf.PdfPages(args.output_name)

    #Separate violinplot for each gene with y-axis over all ontologies
    #y-axis is different ontology
    for gene,g in df.groupby('gene'):
        
        order = g.groupby('annotation')['score'].median().sort_values().index
        
        sns.violinplot(
            x = 'score',
            y = 'annotation',
            orient = 'h',
            order = order,
            data = g,
        )
        plt.ylabel('')
        plt.xlim(-1,1)
        plt.title('{}'.format(gene))
        pdf.savefig(bbox_inches='tight')
        plt.close()

    pdf.close()

if __name__ == '__main__':
    main()

