
bam_ch = Channel.fromPath("/scratch/groups/horence/rob/data/MERFISH_scRNAseq/10X_mapping/merfish_filt_L8TX_1*.bam")

process count_alignments {

    input:
    file(bam) from bam_ch

    output:
    file("num_alignments.txt") into alignment_counts

    script:
    """
    test_SRRS.py \
        --input_bam ${bam}
    """

}

alignment_counts
    .collectFile(
        name: 'all_counts.txt',
        storeDir: 'outputs',
    )

