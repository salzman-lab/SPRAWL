//Process the sample sheet
Channel                                                                         
    .fromPath(params.sample_sheet)                                       
    .splitCsv(header: true)                                                     
    .map{ row ->                                                                
            tuple(                                                              
                row.experiment,                                                     
                row.sample,                                                     
                file(row.path),                                                 
            )                                                                   
    }                                                                           
    .set{ samples_ch } 

//Decide which metrics to score by
metrics_ch = Channel.of('central','peripheral')

process gene_cell_scoring {
    input:
    tuple val(experiment), val(sample), file(hdf5) from samples_ch
    each metric from metrics_ch

    output:
    tuple val(experiment), val(metric), file("${experiment}_${sample}_${metric}.csv") into sample_gene_cell_ch

    script:
    """
    gene_cell_scoring.py \
        --hdf5_path ${hdf5} \
        --metric ${metric} \
        --experiment ${experiment} \
        --sample ${sample} \
        --output_name "${experiment}_${sample}_${metric}.csv" \
        --min_genes_per_cell ${params.min_genes_per_cell} \
        --min_tot_counts_per_cell ${params.min_tot_counts_per_cell} \
    """
}

process merge_gene_cells {
    publishDir "outputs/${params.run_name}/gene_cell", mode: 'copy'

    input:
    tuple val(experiment), val(metric), file('gene_cell') from sample_gene_cell_ch.groupTuple(by: [0,1])

    output:
    tuple val(experiment), val(metric), file("${experiment}_${metric}.csv") into gene_cell_ch

    script:
    """
    awk '(NR == 1) || (FNR > 1)' gene_cell* > "${experiment}_${metric}.csv" 
    """
}


process gene_ont_scoring {
    publishDir "outputs/${params.run_name}/gene_ont", mode: 'copy'

    input:
    tuple val(experiment), val(metric), file(gene_cell) from gene_cell_ch

    output:
    file("${experiment}_${metric}.csv") into gene_ont_ch

    script:
    """
    gene_ont_scoring.py \
        --gene_cell_path ${gene_cell} \
        --output_name "${experiment}_${metric}.csv" \
        --min_cells_per_gene_ont ${params.min_cells_per_gene_ont} \
    """

}


gene_ont_ch.view()



