//Pipeline to segment nuclei from mosaic images and assign them to cells in HDF5 format

params.mosaic_sample_sheet = 'mosaics.csv'
params.fov_sample_sheet = 'fovs.csv'

Channel                                                                         
    .fromPath(params.mosaic_sample_sheet)                                              
    .splitCsv(header: true)                                                     
    .map{ row ->                                                                
            tuple(                                                              
                row.sample,
                row.z,
                file(row.mosaic_path),                                            
                file(row.transform_path),                                            
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


// Nuclei segmentation and assignment to cells
process create_nuclei_boundaries {

    input:
    tuple val(sample), val(z), file(mosaic), file(transform), val(fov), file(hdf5_fov) from mosaic_ch.combine(fov_ch, by: 0)
    
    output:
    file('cell_nuc.gpkg') into nuclei_ch

    script:
    """
    segment_and_assign.py \
        --sample ${sample} \
        --z ${z} \
        --mosaic ${mosaic} \
        --transform ${transform} \
        --fov ${fov} \
        --hdf5_fov ${hdf5_fov} \
    """
}

//Merging all results into a single large geopandas file
//Next step will be to update the SRRS HDF5 files to have nuclei
process merge {
    publishDir 'outputs', mode: 'copy'

    input:
    file('cell_nucs') from nuclei_ch.collect()

    output:
    file('nuclei.gpkg')

    script:
    """
    merge_nuclei.py \
        --cell_nucs cell_nucs* \
    """
}





