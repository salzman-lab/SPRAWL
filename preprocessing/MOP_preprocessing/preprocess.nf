// Inputs downloaded from https://download.brainimagelibrary.org/cf/1c/cf1c1a431ef8d021/processed_data/ July 2021
root_path = "/oak/stanford/groups/horence/rob/isoform_localizations/unprocessed_downloads/dl_scripts/revised_collection"

cell_labels = file("${root_path}/cell_labels.csv")
scriptsDir = "${baseDir}/scripts"
zslices = Channel.of(0, 1, 2, 3, 4, 5, 6)

spot_ch = Channel.fromList([
    //["mouse1sample1",file("${root_path}/spots_mouse1sample1.csv")],
    ["mouse1sample2",file("${root_path}/spots_mouse1sample2.csv")],
    ["mouse1sample3",file("${root_path}/spots_mouse1sample3.csv")],
    ["mouse1sample4",file("${root_path}/spots_mouse1sample4.csv")],
    ["mouse1sample5",file("${root_path}/spots_mouse1sample5.csv")],
    ["mouse1sample6",file("${root_path}/spots_mouse1sample6.csv")],
    ["mouse2sample1",file("${root_path}/spots_mouse2sample1.csv")],
    ["mouse2sample2",file("${root_path}/spots_mouse2sample2.csv")],
    ["mouse2sample3",file("${root_path}/spots_mouse2sample3.csv")],
    ["mouse2sample4",file("${root_path}/spots_mouse2sample4.csv")],
    ["mouse2sample5",file("${root_path}/spots_mouse2sample5.csv")],
    ["mouse2sample6",file("${root_path}/spots_mouse2sample6.csv")],
])

segm_ch = Channel.fromList([
    ["mouse1sample1",file("${root_path}/segmented_cells_mouse1sample1.csv")],
    ["mouse1sample2",file("${root_path}/segmented_cells_mouse1sample2.csv")],
    ["mouse1sample3",file("${root_path}/segmented_cells_mouse1sample3.csv")],
    ["mouse1sample4",file("${root_path}/segmented_cells_mouse1sample4.csv")],
    ["mouse1sample5",file("${root_path}/segmented_cells_mouse1sample5.csv")],
    ["mouse1sample6",file("${root_path}/segmented_cells_mouse1sample6.csv")],
    ["mouse2sample1",file("${root_path}/segmented_cells_mouse2sample1.csv")],
    ["mouse2sample2",file("${root_path}/segmented_cells_mouse2sample2.csv")],
    ["mouse2sample3",file("${root_path}/segmented_cells_mouse2sample3.csv")],
    ["mouse2sample4",file("${root_path}/segmented_cells_mouse2sample4.csv")],
    ["mouse2sample5",file("${root_path}/segmented_cells_mouse2sample5.csv")],
    ["mouse2sample6",file("${root_path}/segmented_cells_mouse2sample6.csv")],
])


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

