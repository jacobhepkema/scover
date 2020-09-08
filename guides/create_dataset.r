#!/usr/bin/env Rscript

suppressMessages(library(AnnotationHub))
suppressMessages(library(assertr))
suppressMessages(library(Biostrings))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Matrix))
suppressMessages(library(optparse))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(scater))
suppressMessages(library(stringr))
options(stringsAsFactors = F)

# Command line argument definitions
option_list <- list(
  make_option(c("-p", "--pool_size"), type="character", default=NULL, 
              help="the amount of cells from the same cell to pool together", 
              metavar="POOL_SIZE"),
  make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="the path to the input SingleCellExperiment RDS file", 
              metavar="DATASET.RDS"),
  make_option(c("-o", "--output"), type="character", default="data.tsv",
              help="the prefix given to the output files (default: 'output_dataset')", 
              metavar="OUTPUT_PREFIX"),
  make_option(c("-c", "--organism"), type="character", default="mouse", 
              help="organism. Options: 'mouse' (default), 'human'", 
              metavar="ORGANISM_NAME"),
  make_option(c("-w", "--output_dir"), type="character", default=".",
              help="output directory. Default: '.'",
              metavar="OUTDIR")
); 

source("scanem_helper_functions.R")

# Parse command line arguments
opt_parser <- optparse::OptionParser(option_list=option_list);
opt <- optparse::parse_args(opt_parser);

# Check if required arguments are supplied
if (is.null(opt$pool_size) |
    is.null(opt$dataset) | 
    is.null(opt$output)){
  optparse::print_help(opt_parser)
  stop("At least these arguments must be supplied:
 1) pool size (-p), 2) dataset (-d)
 3) output prefix (-o)
 Though: keep in mind that the default organism (-c) is 'mouse'\n", call.=FALSE)
}

pool_size <- as.numeric(opt$pool_size)

organism <- base::tolower(opt$organism)
if (!(organism %in% c("mouse", "human"))) {
  stop("Organism should be either 'mouse' or 'human'")
}
cat("Organism:", organism, "\n")

output_prefix <- opt$output
output_directory <- opt$output_dir

assay_name <- base::tolower(opt$assay)
celltype_name <- opt$celltype_name

dataset <- opt$dataset
sce <- readRDS(dataset)

if(!("cell_type1" %in% colnames(colData(sce)))){
  stop("SingleCellExperiment object needs to contain cell type annotation in 
       colData(sce)$cell_type1")
}

# Pool cells using custom function
pool_output <- pool_cells_follow_cells(sce=sce, pool_size = pool_size, 
                                       keep_smaller_pools = FALSE, tries = 1)
pooled_sce <- pool_output[[1]] # Containes pooled sce
pooled_sce_cells_in_pools <- pool_output[[2]] 

if(organism == "mouse"){
  seqs <- Biostrings::readDNAStringSet("./mm10_500up500down_promoters.fa.gz")
} else if(organism == "human") {
  seqs <- Biostrings::readDNAStringSet("./hg38_500up500down_promoters.fa.gz")
} else {
  stop("Invalid organism specified. Options: 'mouse' or 'human'")
}

# Get genes and to-be-predicted data (log-normalized counts of sce)
genelist <- rownames(pooled_sce)
data <- logcounts(pooled_sce) 

# Figure out which genes are not included in the sequences list:
not_in_fasta <- genelist[!(genelist %in% names(seqs))]

genelist <- genelist[genelist %in% names(seqs)]

# Remove cells with over 98% dropouts
cell_dropouts <- colSums(data == 0) / dim(data)[1]
too_high_cell_dropouts <- cell_dropouts > 0.98
data <- data[,!too_high_cell_dropouts]
pooled_sce <- pooled_sce[,!too_high_cell_dropouts]

# Remove genes with over 98% dropout
gene_dropouts <- rowSums(data == 0) / dim(data)[2]
too_high_gene_dropouts <- genelist[gene_dropouts[genelist] > 0.98]
genelist <- genelist[!(gene_dropouts[genelist] > 0.98)]

cat("Genelist of", length(genelist), "genes\n")
cat("Cell/pool number:", ncol(data),"\n")

# Extract sequences
seqs <- base::tolower(as.character(seqs[genelist]))
# Generate comma-separated values to be predicted for each gene
ind <- assertr::col_concat(data[genelist,], sep=",")

dataset <- data.frame(gene=genelist,
                      sequence=seqs,
                      ind=ind)

dataset_name <- paste0(output_directory, "/", output_prefix, "_pool", as.character(pool_size))

write.table(dataset,
            file=gzfile(paste0(dataset_name, ".tsv.gz")),
            sep="\t", row.names = F, col.names = T, quote = F)
write.table(colData(pooled_sce), 
            file=paste0(dataset_name, "_colData.tsv"), 
            sep="\t", row.names = T, col.names = T, quote = F)
write.table(pooled_sce_cells_in_pools, 
            file=paste0(dataset_name, "_cells_in_pools.tsv"), 
            sep="\t", row.names = T, col.names = T, quote = F)

cat("Done\n")