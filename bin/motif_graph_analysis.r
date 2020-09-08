#!/usr/bin/env Rscript

suppressMessages(require(optparse))
suppressMessages(require(stringr))
suppressMessages(require(xtable))
suppressMessages(require(igraph))
options(stringsAsFactors = F)


option_list <- list(
  make_option(c("-d", "--database"), type="character", default=NULL, 
              help="path to the database of motifs used as Tomtom reference", 
              metavar="DATABASE.MEME"), 
  make_option(c("-a", "--tom_alignments"), type="character", default=NULL, 
              help="path to the Tomtom tsv output of motifs aligned to the database", 
              metavar="TOMTOM.TSV"), 
  make_option(c("-m", "--motif_alignments"), type="character", default=NULL,
              help="path to the Tomtom tsv output of motifs aligned to motifs", 
              metavar="TOMTOM_M.TSV"),
  make_option(c("-t", "--threshold"), type="double", default=0.05,
              help="cutoff value of motif matches: significance (q-val)",
              metavar="CUTOFF"),
  make_option(c("-o", "--outdir"), type="character", default=".",
              help="output directory",
              metavar="OUTDIR"))

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- optparse::parse_args(opt_parser);

# Check if required arguments are supplied
if (is.null(opt$database) | 
    is.null(opt$tom_alignments) | 
    is.null(opt$motif_alignments)){
  optparse::print_help(opt_parser)
  stop("At least three arguments must be supplied:
 1) database (-c)\n 2) tom_alignments (-a)\n
 3) motif_alignments (-m)\n", call.=FALSE)
}

# Get full motif names (rather than IDs) from database file
dbfile <- opt$database
db <- readLines(dbfile)
db <- db[str_detect(db, "^MOTIF")]
codes <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][2] })
names <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][3] })
names(names) <- codes

# Get threshold and output directory
threshold <- opt$threshold
outdir <- opt$outdir






# Part one: find motif-alignment clusters and output reproducibility table ========================

# Read motif-motif alignment TSV
motif_alignments <- opt$motif_alignments
all_motif_table <- read.delim(motif_alignments, sep="\t", quote="", header=T,
                              comment.char = "#")
# Remove rows where it aligned to itself
all_motif_table <- all_motif_table[all_motif_table$Query_ID != all_motif_table$Target_ID,]


# Read alignment table and update alignment names to full motif names
tom_alignments <- opt$tom_alignments
all_tom_table <- read.delim(tom_alignments, sep="\t",quote = "", header = T, 
                            comment.char = "#")
all_tom_table$Target_ID <- names[all_tom_table$Target_ID]

# Only use alignments below alignment threshold
all_tom_table <- all_tom_table[all_tom_table$q.value < threshold, ]

# Get edge connections from table
edge_vect <- c()
for(i in 1:nrow(all_tom_table)){
  edge_vect <- c(edge_vect, all_tom_table[i,1], all_tom_table[i,2])
}
edge_weights <- all_tom_table$q.value

# Build graph from motif alignments
g <- graph(edge_vect, directed = F)
bipartite.mapping(g) # Assess bipartite mapping
V(g)$type <- bipartite_mapping(g)$type # Make it a bipartite graph

# Add log-transformed q-values as edge weights
E(g)$weight <- -log2(edge_weights) 

# Set plotting parameters
V(g)$color <- ifelse(V(g)$type, "lightblue", "salmon")
V(g)$shape <- ifelse(V(g)$type, "circle", "square")
E(g)$color <- "lightgray"
V(g)$label.color <- "black" 
V(g)$label.cex <- .8
V(g)$frame.color <-  "gray"
V(g)$size <- 4
E(g)$arrow.size = 0.5
V(g)$size2 = 3

# Plot graph
plot_output_name <- paste0(outdir, "/Motif_alignment_plot.png")
png(plot_output_name, width = 2000, height = 2000)
plot(g, layout = layout_with_graphopt, edge.width=E(g)$weight)
dev.off()

# Cluster graph
g <- simplify(g)
fgc <- igraph::cluster_walktrap(g)
fgc_groups <- igraph::groups(fgc)

# Plot with clusters
plot_output_name <- paste0(outdir, "/Motif_alignment_plot_with_clusters.png")
png(plot_output_name, width = 2000, height = 2000)
plot(fgc, g, layout = layout_with_graphopt, edge.width=E(g)$weight)
dev.off()

# Save clusters to html and tsv
html_file <- paste0(outdir, "/Motif_clusters.html")

n_candidates <- length(unique(sapply(all_tom_table$Query_ID, function(x){return(strsplit(x,split="_")[[1]][1])})))

# Calculate cluster metrics
group_reproducibilities <- sapply(fgc_groups,FUN=function(x, ncan=n_candidates){
  prefixes <- stringr::str_extract(x, "^best|(^[0-9]+)")
  prefixes <- prefixes[!is.na(prefixes)]
  num_exp_found <- length(unique(prefixes))
  return( (num_exp_found/ncan) )
})
cluster_motif_detectors <- sapply(fgc_groups, FUN=function(x, query_options=all_tom_table$Query_ID){
  return(x[x %in% query_options])
})
cluster_not_motif_detectors <- sapply(fgc_groups, FUN=function(x, query_options=all_tom_table$Query_ID){
  return(x[!(x %in% query_options)])
})
cluster_motif_detectors <- sapply(X = cluster_motif_detectors, FUN=function(X){
  return(paste0(X, collapse=" "))
})
cluster_alignments <- sapply(X = cluster_not_motif_detectors, FUN=function(X){
  return(paste0(X, collapse=" "))
})
cluster_df <- data.frame(cluster_motif_detectors, cluster_alignments, cluster_reproducibility=group_reproducibilities)

# Save clusters with metrics to tsv
print(xtable(cluster_df), "html", file = html_file)
tsv_file <- paste0(outdir, "/Motif_clusters.tsv")
write.table(cluster_df, file = tsv_file, sep = "\t", quote = F)

# Save graph object
graph_object_name <- paste0(outdir, "/Motif_alignment_graph_object.RDS")
saveRDS(g, graph_object_name)



cat("\nDone\n")
