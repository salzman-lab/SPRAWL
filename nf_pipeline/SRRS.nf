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
//Default all metrics off, add to list if user turned them on
def metrics_list = []

params.peripheral = false
if(params.peripheral){
    metrics_list.add('peripheral')
}

params.central = false
if(params.central){
    metrics_list.add('central')
}

params.radial = false
if(params.radial){
    metrics_list.add('radial')
}

params.punctate = false
if(params.punctate){
    metrics_list.add('punctate')
}

metrics_ch = Channel.fromList(metrics_list)


//Default params
params.scoring_processes = 5
params.permute_gene_labels = 'no'

process gene_cell_scoring {
    cache 'lenient'
    cpus { params.scoring_processes }

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
        --processes ${params.scoring_processes} \
        --permute_gene_labels ${params.permute_gene_labels} \
    """
}

process merge_gene_cells {
    publishDir "outputs/${params.run_name}/gene_cell", mode: 'copy'

    input:
    tuple val(experiment), val(metric), file('gene_cell') from sample_gene_cell_ch.groupTuple(by: [0,1])

    output:
    tuple val(experiment), val(metric), file("${experiment}_${metric}_gene_cell.csv") into gene_cell_ch

    script:
    """
    awk '(NR == 1) || (FNR > 1)' gene_cell* > "${experiment}_${metric}_gene_cell.csv" 
    """
}

//Duplicate the gene_cell_ch to use for gene_ont and for plotting
gene_cell_ch.into{ scoring_gene_cell_ch; plotting_gene_cell_ch }


process gene_ont_scoring {
    publishDir "outputs/${params.run_name}/gene_ont", mode: 'copy'

    input:
    tuple val(experiment), val(metric), file(gene_cell) from scoring_gene_cell_ch

    output:
    file("${experiment}_${metric}_gene_ont.csv") into gene_ont_ch

    script:
    """
    gene_ont_scoring.py \
        --gene_cell_path ${gene_cell} \
        --output_name "${experiment}_${metric}_gene_ont.csv" \
        --min_cells_per_gene_ont ${params.min_cells_per_gene_ont} \
    """

}



process gene_cell_plotting {
    publishDir "outputs/${params.run_name}/plots", mode: 'copy'

    input:
    tuple val(experiment), val(metric), file(gene_cell) from plotting_gene_cell_ch

    output:
    file("${experiment}_${metric}_gene_cell_plots.pdf") into gene_cell_plots_ch

    script:
    """
    gene_cell_plotting.py \
        --gene_cell_path ${gene_cell} \
        --metric ${metric} \
        --experiment ${experiment} \
        --output_name "${experiment}_${metric}_gene_cell_plots.pdf" \
    """
}



gene_ont_ch.view()

