// Inputs downloaded from https://download.brainimagelibrary.org/cf/1c/cf1c1a431ef8d021/processed_data/ July 2021

scriptsDir = "${baseDir}/scripts"

// User should set these items
sample_name = "Liver"
spot_path = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/barcodes.csv"
boundary_paths = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/features/*_11*.hdf5"

// Start
spot_ch = Channel.fromPath(spot_path)
boundary_ch = Channel.fromPath(boundary_paths)


// Restructure the points into an RTree index 
process create_spot_rtree {

    input:
    file(spot_path) from spot_ch

    output:
    file("rtree.idx") into rtree_idx_ch
    file("rtree.dat") into rtree_dat_ch

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

    input:
    file(boundary_hdf5) from boundary_ch
    file(idx) from rtree_idx_ch
    file(dat) from rtree_dat_ch

    output:
    file('spot_assigned_fov.hdf5') into hdf5_fovs

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


hdf5_fovs.view()


/*
// Merge fovs
process merge_fovs {
    publishDir 'finished_outputs', mode: 'copy'

    input:
    tuple val(mouse_sample),file('hdf5_slice*') from zslice_cells_ch.groupTuple()

    output:
    tuple val(mouse_sample),file("${mouse_sample}.hdf5") into final_out_ch

    script:
    """                                                                         
    ${scriptsDir}/merge_slices.py \
        --slices hdf5_slice* \
        --cell_labels ${cell_labels} \
        --mouse_sample ${mouse_sample} \
    """  
}


final_out_ch.view()
*/
