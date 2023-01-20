library(Gviz)
library(GenomicFeatures)
library(data.table)

#support for plus/minus strand bam plotting from:
#https://support.bioconductor.org/p/56047/
library(Rsamtools)
strandedBamImport <- function (file, selection) {
    if (!file.exists(paste(file, "bai", sep = ".")))
        stop("Unable to find index for BAM file '", file, "'. You can
build an index using the following command:\n\t",
             "library(Rsamtools)\n\tindexBam(\"", file, "\")")
    sinfo <- scanBamHeader(file)[[1]]
    res <- if (!as.character(seqnames(selection)[1]) %in%
               names(sinfo$targets)) {
        mcols(selection) <- DataFrame(score = 0)
        selection
    }else {
        param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                              which = selection, flag =
scanBamFlag(isUnmappedQuery = FALSE))
        x <- scanBam(file, param = param)[[1]]
        gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]],
width = x[["qwidth"]]), seqnames=seqnames(selection)[1])
        grs <- split(gr, strand(gr))
        cov <- lapply(grs[c("+", "-")], function(y)
coverage(ranges(y),
width=end(selection)))
        pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y),
end(y))))))
        if(length(pos)==0){
            mcols(selection) <- DataFrame(plus=0, minus=0)
            selection
        }else{
            GRanges(seqnames = seqnames(selection)[1],
ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)),
                    plus=as.numeric(cov[["+"]][head(pos, -1)]),
                    minus=-as.numeric(cov[["-"]][head(pos, -1)]))
        }
    }
    return(res)
}
                                         

args = commandArgs(trailingOnly=TRUE)

gene_name = args[1]

txdb_path = "txdb.mm10.sqlite"
bed_path = "MERFISH_genes.bed"
srrs_path = "/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/outputs/gene_ontology/MOp_periphal_ReadZs_10X_gene_ontology.csv"
bam_stem = "/scratch/groups/horence/rob/data/MERFISH_scRNAseq/10X_mapping/merged_by_celltype/"
                                         
bam_paths = Sys.glob(paste(bam_stem,"*.bam",sep=""))
bam_ontologies = lapply(lapply(bam_paths,basename),tools::file_path_sans_ext)
                                         
#find which cell-types to plot from SRRS
srrs_df = fread(srrs_path)
srrs_df = srrs_df[gene == gene_name] #filter to results from the gene of interest
srrs_df = srrs_df[bh_p < 0.05] #filter by bh_p
srrs_df = srrs_df[order(med_score)] #order by effect size
srrs_df = srrs_df[ontology %in% bam_ontologies] #keep only ontologies which have corresponding bams
srrs_df = unique(srrs_df,by='ontology') #deduplicate by ontology

n_onts_to_plot = 6
                                         
if(nrow(srrs_df) == 0){
    print("No matching cell-types")
    stop()
} else if (nrow(srrs_df) <= n_onts_to_plot){
    #plot all cell types
    n_onts_to_plot = nrow(srrs_df)
    cell_types = srrs_df$ontology
    effect_sizes = srrs_df$med_score
} else {
    #plot top n/2 and bottom n/2 cell-types by med score (effect size)
    n = n_onts_to_plot/2
    cell_types = c(head(srrs_df,n)$ontology, tail(srrs_df,n)$ontology)
    effect_sizes = c(head(srrs_df,n)$med_score, tail(srrs_df,n)$med_score)
}

                                         
gene_locs = fread(bed_path)
gene_row = gene_locs[gene == gene_name]

chr = gene_row[, get("#chr")]
window_start = gene_row$start
window_end   = gene_row$end
outfile = paste(gene_name,".png",sep="")
tracks = list()
                                         
####################################
#                                  #
#                                  #
#            Gene Track            #
#                                  #
#                                  #
####################################
txdb = loadDb(txdb_path)

geneTrack = GeneRegionTrack(
    range = txdb,
    chromosome = chr,
    start = window_start,
    end = window_end,
    name = gene_name,
)

geneTrack <- geneTrack[
    feature(geneTrack) %in% c("CDS", "UTR")
]
tracks = append(tracks,geneTrack)

gTrack <- GenomeAxisTrack()                                         
tracks = append(tracks,gTrack)

####################################
#                                  #
#                                  #
#            Data Track            #
#                                  #
#                                  #
####################################
for(i in 1:n_onts_to_plot){
    ct = cell_types[[i]]
    ct = sub(" ","_",ct)
    
    effect_size = effect_sizes[[i]]
    effect_size = format(round(effect_size, 2), nsmall = 2)
    
    cov_track = DataTrack(
        range = paste(bam_stem,ct,".bam",sep=""),
        name = paste(ct,effect_size,sep="\n"),
        chromosome = chr,
        importFunction=strandedBamImport,
        stream=TRUE,
    )
    tracks = append(tracks,cov_track)
}

####################################
#                                  #
#                                  #
#            Plotting              #
#                                  #
#                                  #
####################################
png(outfile, width = 800, height = 800)
par(mar=c(5,6,4,1)+.1)
plotTracks(
    tracks,
    type = "h",
    from = window_start, 
    to = window_end,
    thinBoxFeature = "UTR",
    collapse = F,
    transcriptAnnotation = "gene",
    col=c("red","blue"),
    groups=c("+", "-"),
    fontcolor.exon = "black",
    legend = F,
    cex = 0.6
)

dev.off()
