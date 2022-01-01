library(Gviz)
library(GenomicFeatures)

args = commandArgs(trailingOnly=TRUE)

print(args)

chr = args[1]
window_start = strtoi(args[2])
window_end   = strtoi(args[3])
outfile = args[4]
txdb_path = args[5]
celltypes = args[6]


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
)

geneTrack <- geneTrack[
    feature(geneTrack) %in% c("CDS", "UTR")
]

gTrack <- GenomeAxisTrack()

####################################
#                                  #
#                                  #
#            Data Tracks           #
#                                  #
#                                  #
####################################
data_tracks = c()

for (i in 5:length(args)){
    print(args[i])
    data_track <-DataTrack(
        range = args[i],
        name=basename(args[i]),
        chromosome=chr,
    )
    data_tracks <- c(data_tracks, data_track)
}


####################################
#                                  #
#                                  #
#            Plotting              #
#                                  #
#                                  #
####################################
png(outfile, width = 800, height = 400)

all_tracks = c(gTrack, data_tracks, geneTrack)
    

plotTracks(
    all_tracks,
    type = "h",
    from = window_start, 
    to = window_end,
    thinBoxFeature = "UTR",
    collapse = F,
    transcriptAnnotation = "gene",
    fontcolor.exon = "black",
    legend = T,
    cex = 0.6
)

dev.off()

