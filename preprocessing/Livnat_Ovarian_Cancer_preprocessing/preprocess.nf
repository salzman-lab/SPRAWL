
scriptsDir = "${baseDir}/scripts"

spots_ch = Channel.fromFilePairs(
    "/scratch/users/rbierman/spatial_collaboration/detected_transcripts/*_detected_transcripts.csv",
    size: 1,
    flat: true,
)

boundaries = Channel.fromFilePairs(
    "/scratch/users/rbierman/spatial_collaboration/cell_boundaries/*.hdf5",
    size: -1,
){ file -> file.name.replaceAll(/_feature_data_.*$/,'') }

min_gene_spots = 2
min_tot_spots = 100
min_uniq_genes = 10


// Restructure the points into an RTree index 
process create_spot_rtree {

    memory '20 GB'
    time '48h'

    input:
    tuple val(sample), file(spot_path) from spots_ch

    output:
    tuple val(sample), file("rtree.idx"),file("rtree.dat") into rtree_ch

    script:
    """
    ${scriptsDir}/index_spots.py \
        --spot_path ${spot_path} \
    """

    stub:
    """
    echo ${spot_path}
    touch rtree.idx
    touch rtree.dat
    """
}

// Assign spots to cells
process assign_spots_to_cells {
    cache 'lenient'

    //Input is a bit ridiculous.
    //Basically merging boundary and rtrees by sample name and flattening them
    //must be a better way to do this
    //vals sample1 and sample2 are the same, but I don't know how to avoid duplication
    input:
    tuple val(sample1),file("rtree.idx"),file("rtree.dat"),val(sample2),file(boundary_hdf5) from rtree_ch.cross(boundaries.transpose()).map{t -> t.flatten()}

    output:
    tuple val(sample1),file('spot_assigned_fov.hdf5') into hdf5_fovs

    script:
    """
    ${scriptsDir}/assign_spots_to_cells.py \
        --hdf5_path ${boundary_hdf5} \
        --rtree_name rtree \
    """

    stub:
    """
    echo ${boundary_hdf5}
    touch spot_assigned_fov.hdf5
    """
}


// Filter and calculate gene-cell variances
process filter_and_variance_calc {
    cache 'lenient'
    memory '20 GB'
    time '1h'
    cpus 10 


    input:
    tuple val(sample1),file(srrs_hdf5) from hdf5_fovs

    output:
    tuple val(sample1),file('spot_assigned_fov_gv.hdf5') into hdf5_gv_fovs

    script:
    """
    ${scriptsDir}/filter_gene_cell_variance.py \
        --hdf5_path ${srrs_hdf5} \
        --out_path spot_assigned_fov_gv.hdf5 \
        --min_gene_spots ${min_gene_spots} \
        --min_tot_spots ${min_tot_spots} \
        --min_genes ${min_uniq_genes} \
    """

    stub:
    """
    touch spot_assigned_fov_gv.hdf5
    """
}


// Merge fovs
process merge_fovs {
    publishDir 'finished_outputs', mode: 'copy'

    //Kind of like collectFile()
    input:
    tuple val(sample), file('hdf5_fov') from hdf5_gv_fovs.groupTuple()

    output:
    file("${sample}.hdf5") into final_out_ch

    script:
    """                                                                         
    ${scriptsDir}/merge_fovs.py \
        --hdf5_fovs hdf5_fov* \
        --sample ${sample} \
    """  

    stub:
    """
    touch ${sample}.hdf5
    """

}


final_out_ch.view()

