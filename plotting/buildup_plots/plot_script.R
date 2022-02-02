#Automating gene plotting
#Need a table of gene/ontology with mean RZS score

library(Gviz)
library(GenomicFeatures)
library(data.table)

options(ucscChromosomeNames=FALSE)

args = commandArgs(trailingOnly=TRUE)
gene_name = args[1]
gene_padding = 0.05 #how much to pad upstream and downstream of the gene

txdb_path = "txdb.mm10.sqlite"
bed_path = "MERFISH_genes.bed"
bam_stem = "/scratch/groups/horence/rob/data/MERFISH_scRNAseq/10X_mapping/merged_by_celltype/"
rzs_score_path = "/oak/stanford/groups/horence/rob/isoform_localizations/SRRS/analysis/gene_ontology_rzs.csv"
bam_paths = Sys.glob(paste(bam_stem,"*.bam",sep=""))
                                         
gene_locs = fread(bed_path)
gene_row = gene_locs[gene == gene_name]

gene_buffer = (gene_row$end-gene_row$start)*gene_padding

chr = gene_row[, get("#chr")]
window_start = round(gene_row$start-gene_buffer)
window_end   = round(gene_row$end+gene_buffer)
gene_strand  = gene_row$strand
outfile = paste(gene_name,".png",sep="")
tracks = list()

df = fread(rzs_score_path)
df = df[gene == gene_name]
df = df[order(mean_rzs)]

cell_types = df$ontology

stopifnot(length(cell_types) > 0)
 
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
    strand = gene_strand,
    name = gene_name,
)

geneTrack <- geneTrack[
    gene(geneTrack) == gene_name
]
tracks = append(tracks,geneTrack)

####################################
#                                  #
#                                  #
#            Data Tracks           #
#                                  #
#                                  #
####################################
data_tracks = list()

for(ct in cell_types){
        
    cov_track = AlignmentsTrack(
        range = paste(bam_stem,ct,".bam",sep=""),
        name = ct,
        chromosome = chr,
        #col = colors[[i]],
        #fill = colors[[i]],
        strand = gene_strand,
        #legend = TRUE,
        type = "coverage",
        window = -1,
        windowSize = 1,
    )
    data_tracks = append(data_tracks,cov_track)
}
tracks = append(tracks, data_tracks)

####################################
#                                  #
#                                  #
#            Plotting              #
#                                  #
#                                  #
####################################
png(outfile, width = 800, height = 400)
par(mar=c(5,6,4,1)+.1)

plotTracks(
    tracks,
    sizes = rep(1/length(tracks),length(tracks)), #give equal vertical space to each plot
    from = window_start, 
    to = window_end,
    thinBoxFeature = c("utr3","utr5"),
    collapse = F,
    transcriptAnnotation = "gene",
    fontcolor.exon = "black"
)

dev.off()
