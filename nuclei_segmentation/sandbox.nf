//Pipeline to segment nuclei from mosaic images and assign them to cells in HDF5 format

scriptsDir = "${baseDir}/scripts"

params.mosaic_sample_sheet = 'mosaics.csv'
params.fov_sample_sheet = 'fovs_subset.csv'

Channel                                                                         
    .fromPath(params.mosaic_sample_sheet)                                              
    .splitCsv(header: true)                                                     
    .map{ row ->                                                                
            tuple(                                                              
                row.sample,
                row.z,
                file(row.path),                                            
            )                                                                   
    }                                                                           
    .set{ mosaic_ch }  

Channel                                                                         
    .fromPath(params.fov_sample_sheet)                                              
    .splitCsv(header: true)                                                     
    .map{ row ->                                                                
            tuple(                                                              
                row.sample,
                row.fov,
                file(row.path),                                            
            )                                                                   
    }                                                                           
    .set{ fov_ch }  



// Nuclei segmentation and assignment
process create_nuclei_boundaries {

    input:
    tuple val(sample), val(z), file(mosaic), val(fov), file(hdf5_fov) from mosaic_ch.combine(fov_ch, by: 0)
    
    output:
    tuple val(sample), val(z), val(fov), file("cells_with_nuclei.hdf5") into nuclei_ch

    script:
    """
    """

    stub:
    """
    touch nuclei.gpkg
    """
}

// Assign nuclei to cells
process assign_nuclei_to_cells {

    input:
    tuple val(sample), val(z), val(fov), file("nuclei.gpkg") from nuclei_ch

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


