// Inputs downloaded from https://download.brainimagelibrary.org/cf/1c/cf1c1a431ef8d021/processed_data/ July 2021

scriptsDir = "${baseDir}/scripts"

// User should set these items

// LIVER
//sample = "liver"
//spot_path = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/barcodes.csv"
//boundary_paths = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/features/*.hdf5"

// KIDNEY 1
//sample = "kidney_111921"
//spot_path = "data/VZG116/MsKidney_CellBoundary_VZG116_111921/barcodes.csv"
//boundary_paths = "data/VZG116/MsKidney_CellBoundary_VZG116_111921/features/*.hdf5"

// KIDNEY 2
sample = "kidney_121021"
spot_path = "data/VZG116/MsKidney_CellBoundary_VZG116_121021/barcodes.csv"
boundary_paths = "data/VZG116/MsKidney_CellBoundary_VZG116_121021/features/*.hdf5"


// Start
spot_ch = Channel.fromPath(spot_path)
boundary_ch = Channel.fromPath(boundary_paths)


// Restructure the points into an RTree index 
process create_spot_rtree {

    input:
    file(spot_path) from spot_ch

    output:
    tuple file("rtree.idx"),file("rtree.dat") into rtree_idx_ch

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

    input:
    tuple file("rtree.idx"),file("rtree.dat") from rtree_idx_ch
    each file(boundary_hdf5) from boundary_ch

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


// Merge fovs
process merge_fovs {
    publishDir 'finished_outputs', mode: 'copy'

    input:
    file('hdf5_fov') from hdf5_fovs.collect()

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
    echo ${boundary_hdf5}
    touch ${sample}.hdf5
    """

}


final_out_ch.view()
