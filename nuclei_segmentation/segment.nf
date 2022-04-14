//Pipeline to segment nuclei from mosaic images and assign them to cells in HDF5 format

scriptsDir = "${baseDir}/scripts"

min_filt_spots = 100
min_uniq_genes = 10

// TODO expand to process all LiverXSliceY pairs
mosaic_ch = Channel.fromFilePairs("/scratch/groups/horence/rob/data/vz_liver_showcase/Liver1Slice1/images/mosaic_DAPI_z*.tif")
fov_ch = Channel.fromFilePairs("/scratch/groups/horence/rob/data/vz_liver_showcase/Liver1Slice1/cell_boundaries/feature_data_*.hdf5")

mosaic_ch.view()

// Nuclei segmentation
process create_nuclei_boundaries {

    input:
    tuple val(z), file(mosaic) from mosaic_ch
    tuple val(fov),file(hdf5_fov) from fov_ch
    

    output:
    tuple val(z), val(fov), file("nuclei.gpkg") into nuclei_ch

    script:
    """
    """

    stub:
    """
    touch nuclei.gpkg
    """
}

// Assign spots to cells
process assign_nuclei_to_cells {
    cache 'lenient'

    input:
    tuple val(z), val(fov), file("nuclei.gpkg") from nuclei_ch

    output:
    tuple val(fov),file('nuclei_assigned_fov.hdf5') into nuclei_hdf5_ch

    script:
    """
    """

    stub:
    """
    touch nuclei_assigned_fov.hdf5
    """
}


// Merge fovs for each Sample/Slice
process merge_fovs {
    publishDir 'finished_outputs', mode: 'copy'

    //Kind of like collectFile()
    input:
    tuple val(fov), file('hdf5_fov') from nuclei_hdf5_ch

    output:
    file("temp.hdf5") into final_out_ch

    script:
    """                                                                         
    """  

    stub:
    """
    touch temp.hdf5
    """

}


final_out_ch.view()

