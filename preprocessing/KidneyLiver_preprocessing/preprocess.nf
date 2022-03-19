// Inputs downloaded from https://download.brainimagelibrary.org/cf/1c/cf1c1a431ef8d021/processed_data/ July 2021

scriptsDir = "${baseDir}/scripts"

stem_path = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021"

spots_path = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/barcodes.csv"
boundary_paths = "data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/features/*.hdf5"

spot_ch = Channel.fromList([
    ["MsKidney_1",file("data/VZG116/MsKidney_CellBoundary_VZG116_111921/barcodes.csv")],
    ["MsKidney_2",file("data/VZG116/MsKidney_CellBoundary_VZG116_121021/barcodes.csv")],
    ["MsLiver",file("data/VZG116/MsLiver_Cellbound_VZG116_V1_JH_09-18-2021/barcodes.csv")],
])

segm_ch = Channel.fromPath(boundary_paths)


// Restructure the points into an RTree index 
process create_spot_rtree {

    input:
    tuple val(mouse_sample),file(spot_path) from spot_ch
    each zslice from zslices

    output:
    tuple val(mouse_sample),val(zslice),file('rtree.idx'),file('rtree.dat') into rtree_ch

    script:
    """
    ${scriptsDir}/index_spots.py \
        --spot_path ${spot_path} \
        --z ${zslice} \
    """
}

// Cross the RTree outpus with the segment files
// The following code results in a channel with the following tuples:
//      Sample              Segment path               Z    IDX         DAT
//  ---------------------------------------------------------------------------
//  [mouse1sample1, segmented_cells_mouse1sample1.csv, 1, rtree.idx, rtree.dat],
//  [mouse1sample1, segmented_cells_mouse1sample1.csv, 0, rtree.idx, rtree.dat],
segm_ch
    .cross(rtree_ch)
    .flatMap{ t1,t2 -> collect{ [t1[0],t1[1],t2[1],t2[2],t2[3]] }}
    .set{ segm_spots }


// Assign spots to cells
process assign_spots_to_cells {

    input:
    tuple val(mouse_sample),file(segm_path),val(zslice),file(idx),file(dat) from segm_spots

    output:
    tuple val(mouse_sample),file('slice_cells.hdf5') into zslice_cells_ch

    script:
    """
    ${scriptsDir}/assign_spots_to_cells.py \
        --segm_path ${segm_path} \
        --rtree_name rtree \
        --z ${zslice} \
    """
}


// Merge slices
process merge_slices {
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

