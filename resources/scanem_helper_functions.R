get_sce_sparsity <- function(sce){
  # Computes sparsity of the counts matrix of a SingleCellExperiment
  require(Matrix)
  # (double divide rather than divided by (row*col) so it can handle large 
  #  matrices without integer overflow)
  return(1 - (Matrix::nnzero(counts(sce))/nrow(sce)/ncol(sce)))
}

get_dataset_sparsity <- function(tsv_file){
  # Computes sparsity of the to-be-predicted values of a given scanem dataset
  require(stringr)
  
  # read dataset and initialize array
  dataset <- read.csv(tsv_file, header=T, sep="\t", quote="")
  cells <- length(stringr::str_split(dataset$ind[1], ",")[[1]])
  vals_array <- matrix(nrow = nrow(dataset), ncol=cells)
  
  # cycle through data points, store in array
  for(i in 1:nrow(vals_array)){
    vals_array[i,] <- str_split(dataset$ind[i], ",")[[1]]
  }
  
  # calculate dataset sparsity
  dataset_sparsity <- sum(vals_array == 0)/nrow(vals_array)/ncol(vals_array)
  # (double divide rather than divided by (row*col) so it can handle large 
  #  matrices without integer overflow)
  
  return(dataset_sparsity)
}

shan_ent <- function(vect){
  vect <- abs(vect)
  ent_total <- 0
  for(i in 1:length(vect)){
    if(!(vect[i] > 0)){   # if vect_i = 0: set to 0
      ent_total <- ent_total + 0
    } else {
      ent_total <- ent_total + (vect[i]/sum(vect) * log2(vect[i]/sum(vect)))
    }
  }
  return((-1*ent_total))
}

# function for pooling
pool_cells <- function(sce, pool_size, tries=1, keep_smaller_pools=F, max_pools=NULL){
  # cell pooling 
  # requires single cell experiment, pool size
  # note: sce cell types have to be stored in sce$cell_type1
  # optional: tries. will get try with lowest dropout ratio
  
  require(SingleCellExperiment)
  require(scater)
  require(stringr)
  
  options(stringsAsFactors = F)
  
  cell_types <- as.character(unique(sce$cell_type1[order(sce$cell_type1)]))
  sce <- sce[,order(sce$cell_type1)]
  
  use_max_pools <- FALSE
  if(!is.null(max_pools)){
    cat("Using max pools / cell type", max_pools, "\n")
    use_max_pools <- TRUE
  }
  
  cat("Pooling cells per cell type with pool_size", pool_size, "\n")
  
  best_dropout_ratio <- 1
  for(try in 1:tries){
    if(tries > 1) { cat("Try", try, "\n") }
    total_new_counts <- cbind()
    
    for(i in 1:length(cell_types)){
      curr_cell_type <- cell_types[i]
      curr_data <- sce[,sce$cell_type1 == curr_cell_type]
      
      curr_cells <- ncol(curr_data)
      
      include_celltype <- TRUE
      
      if(!keep_smaller_pools){
        if(curr_cells < pool_size){
          include_celltype <- FALSE
        } else {
          new_ncols <- curr_cells - (curr_cells %% pool_size)
          curr_data <- curr_data[,sample(new_ncols)]
        }
      }
      
      if(include_celltype){
        # random column shuffle:
        curr_data <- curr_data[,sample(ncol(curr_data))]
        curr_counts <- counts(curr_data)
        
        # sum per pool_size columns
        curr_fragmented_counts <- sapply(seq(1,ncol(curr_counts), 
                                             by=pool_size), function(z) {
                                               indx <- z:(z+(pool_size-1))
                                               rowSums(curr_counts[,indx[indx <= ncol(curr_counts)], drop=F])})
        
        colnames(curr_fragmented_counts) <- paste0(paste0(curr_cell_type, "_pool"), 
                                                   1:ncol(curr_fragmented_counts))
        
        if(use_max_pools){
          num_cols <- ifelse(ncol(curr_fragmented_counts) > max_pools, 
                             max_pools, 
                             ncol(curr_fragmented_counts))
          curr_fragmented_counts <- curr_fragmented_counts[,c(1:num_cols),drop=F]
        }
        
        total_new_counts <- cbind(total_new_counts, curr_fragmented_counts)
      } else {
        cat("Not including cell type", curr_cell_type, "\n")
      }
    }
    
    celltypes_combined <- stringr::str_remove(colnames(total_new_counts), pattern = "_pool[0-9]+$")
    combined_sce <- SingleCellExperiment(assays=list(counts=total_new_counts), 
                                         colData=data.frame(cell_type1=celltypes_combined),
                                         rowData=rowData(sce))
    
    combined_sce <- scater::logNormCounts(combined_sce)
    
    dropout_ratio <- get_sce_sparsity(combined_sce)
    cat("\nDropout ratio:", round(dropout_ratio, 3), "\n")
    if(dropout_ratio < best_dropout_ratio){
      best_dropout_ratio <- dropout_ratio
      best_sce <- combined_sce
    }
  }
  
  cat("Cell number:", ncol(combined_sce), "\n")
  cat("Done\n\n")
  
  return(best_sce)
}

get_cluster_names <- function(cluster_df){
  cluster_names <- sapply(cluster_df$cluster_alignments, FUN=function(x){
    if(str_detect(x, "\\(")){
      if(length(str_split(x, " ")[[1]]) > 1){
        if(length(str_split(x, " ")[[1]]) == sum(str_detect(str_split(x, " ")[[1]], "\\("))){
          curr_names <- str_remove_all(sapply(str_split(str_split(
            cluster_df$cluster_alignments[6], " ")[[1]], "_"), FUN=function(x) { return(x[1]) }), "\\(|\\)")
          if(length(curr_names) <= 5){
            names <- curr_names
          } else {
            names <- c(curr_names[1:5], "...")
          }
        } else {
          curr_names <- str_split(x, " ")[[1]]
          curr_names <- curr_names[!str_detect(curr_names, "\\(")]
          if(length(curr_names) <= 5){
            names <- curr_names
          } else {
            names <- c(curr_names[1:5], "...")
          }
        }
      } else {
        curr_names <- str_split(x, " ")[[1]]
        names <- curr_names
      }
    } else {
      curr_names <- str_split(x, " ")[[1]]
      if(length(curr_names) <= 5){
        names <- curr_names
      } else {
        names <- c(curr_names[1:5], "...")
      }
    }
    print(paste(names, collapse = "/"))
    return(paste(names, collapse = "/"))
  })
  return(cluster_names)
}


pool_cells_follow_cells <- function(sce, pool_size, tries=1, keep_smaller_pools=F){
  # cell pooling 
  # requires single cell experiment, pool size
  # note: sce cell types have to be stored in sce$cell_type1
  # optional: tries. will get try with lowest dropout ratio
  
  require(SingleCellExperiment)
  require(scater)
  require(stringr)
  
  options(stringsAsFactors = F)
  
  cell_types <- as.character(unique(sce$cell_type1[order(sce$cell_type1)]))
  sce <- sce[,order(sce$cell_type1)]
  
  cat("Pooling cells per cell type with pool_size", pool_size, "\n")
  
  best_dropout_ratio <- 1
  for(try in 1:tries){
    if(tries > 1) { cat("Try", try, "\n") }
    total_new_counts <- cbind()
    total_new_barcodes <- rbind()
    
    for(i in 1:length(cell_types)){
      curr_cell_type <- cell_types[i]
      curr_data <- sce[,sce$cell_type1 == curr_cell_type]
      
      curr_barcodes <- colnames(curr_data)
      curr_cells <- ncol(curr_data)
      
      include_celltype <- TRUE
      
      if(!keep_smaller_pools){
        if(curr_cells < pool_size){
          include_celltype <- FALSE
        } else {
          new_ncols <- curr_cells - (curr_cells %% pool_size)
          curr_sample <- sample(new_ncols)
          curr_data <- curr_data[,curr_sample]
          curr_barcodes <- curr_barcodes[curr_sample]
        }
      }
      
      if(include_celltype){
        # random column shuffle:
        curr_sample <- sample(ncol(curr_data))
        
        curr_data <- curr_data[,curr_sample]
        curr_counts <- counts(curr_data)
        curr_barcodes <- curr_barcodes[curr_sample]
        curr_fragmented_barcodes <- sapply(seq(1,ncol(curr_counts), 
                                               by=pool_size), function(z) {
                                                 return(curr_barcodes[z:(z+(pool_size-1))])
                                               })
        
        curr_fragmented_barcodes <- t(curr_fragmented_barcodes)
        
        # sum per pool_size columns
        curr_fragmented_counts <- sapply(seq(1,ncol(curr_counts), 
                                             by=pool_size), function(z) {
                                               indx <- z:(z+(pool_size-1))
                                               rowSums(curr_counts[,indx[indx <= 
                                                                           ncol(curr_counts)], 
                                                                   drop=F])})
        
        colnames(curr_fragmented_counts) <- paste0(paste0(curr_cell_type, "_pool"), 
                                                   1:ncol(curr_fragmented_counts))
        
        rownames(curr_fragmented_barcodes) <- paste0(paste0(curr_cell_type, 
                                                            "_pool"), 
                                                     1:ncol(curr_fragmented_counts))
        
        total_new_counts <- cbind(total_new_counts, 
                                  curr_fragmented_counts)
        total_new_barcodes <- rbind(total_new_barcodes, 
                                    curr_fragmented_barcodes)
      } else {
        cat("Not including cell type", curr_cell_type, "\n")
      }
    }
    
    celltypes_combined <- stringr::str_remove(colnames(total_new_counts), 
                                              pattern = "_pool[0-9]+$")
    combined_sce <- SingleCellExperiment(assays=list(counts=total_new_counts), 
                                         colData=data.frame(
                                           cell_type1=celltypes_combined),
                                         rowData=rowData(sce))
    
    combined_sce <- scater::logNormCounts(combined_sce)
    
    dropout_ratio <- get_sce_sparsity(combined_sce)
    cat("\nDropout ratio:", round(dropout_ratio, 3), "\n")
    if(dropout_ratio < best_dropout_ratio){
      best_dropout_ratio <- dropout_ratio
      best_sce <- combined_sce
    }
  }
  
  cat("Cell number:", ncol(best_sce), "\n")
  cat("Done\n\n")
  
  return(list(best_sce, total_new_barcodes))
}

to_z <- function(df){
  rownames_df <- rownames(df)
  colnames_df <- colnames(df)
  df_mat <- as.matrix(df)
  for(i in 1:nrow(df_mat)){
    df_mat[i,] <- (df_mat[i,] - mean(df_mat[i,]))/sd(df_mat[i,])
  }
  colnames(df_mat) <- colnames_df
  rownames(df_mat) <- rownames_df
  return(as.data.frame(df_mat))
}

get_max_motif_activations <- function(M){
  # rectify motifs:
  M <- apply(M, c(1,2,3), FUN=function(x){ ifelse(x>0,x,0) })
  
  d <- dim(M)[3] # amount of motifs
  m <- dim(M)[2] # length of motif
  
  motif_sums <- c()
  pb <- txtProgressBar(0,d,style=3)
  for(i in 1:d){
    setTxtProgressBar(pb,i)
    motif_sum <- 0
    for(j in 1:m){
      if(sum(M[,j,i]) != 0){
        motif_sum <- motif_sum + M[which(M[,j,i] == max(M[,j,i])),j,i]
      }
    }
    motif_sums <- c(motif_sums, motif_sum)
  }
  
  return(motif_sums)
}


pool_sce_given_groups <- function(sce, groups_file){
  require(SingleCellExperiment)
  require(scater)
  require(stringr)
  
  options(stringsAsFactors = F)
  
  get_sce_sparsity <- function(sce){
    # Computes sparsity of the counts matrix of a SingleCellExperiment
    require(Matrix)
    # (double divide rather than divided by (row*col) so it can handle large 
    #  matrices without integer overflow)
    return(1 - (Matrix::nnzero(counts(sce))/nrow(sce)/ncol(sce)))
  }
  
  cell_types <- as.character(unique(sce$cell_type1[order(sce$cell_type1)]))
  sce <- sce[,order(sce$cell_type1)]
  
  groups_file <- t(read.csv(groups_file, sep="\t", header=F, row.names = 1))
  rownames(groups_file) <- 1:nrow(groups_file)  
  total_new_counts <- cbind()
  for(i in 1:ncol(groups_file)){
    curr_pool_name <- colnames(groups_file)[i]
    cells_in_pool <- groups_file[,curr_pool_name]
    curr_counts <- counts(sce[,cells_in_pool])
    total_new_counts <- cbind(total_new_counts, rowSums(curr_counts))
  }
  
  colnames(total_new_counts) <- colnames(groups_file)
  celltypes_combined <- stringr::str_remove(colnames(total_new_counts), 
                                            pattern = "_pool[0-9]+$")
  
  combined_sce <- SingleCellExperiment(assays=list(counts=total_new_counts), 
                                       colData=data.frame(cell_type1=
                                                            celltypes_combined),
                                       rowData=rowData(sce))
  
  combined_sce <- scater::logNormCounts(combined_sce)
  
  dropout_ratio <- get_sce_sparsity(combined_sce)
  cat("\nDropout ratio:", round(dropout_ratio, 3), "\n")
  cat("Cell number:", ncol(combined_sce), "\n")
  cat("Done\n\n")
  
  return(combined_sce)
}

get_tomtom_names <- function(dbfile){
  require(stringr)
  
  db <- readLines(dbfile)
  db <- db[str_detect(db, "^MOTIF")]
  codes <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][2] })
  tomtom_names <- sapply(db, FUN=function(x){ str_split(x, " ")[[1]][3] })
  names(tomtom_names) <- codes
  return(tomtom_names)
}

read_tomtom_output <- function(motif_db_path, exp_folder, 
                               alignment_threshold = 0.05){
  tomtom_path <- paste0(exp_folder, "/all_tomtom/tomtom.tsv")
  
  # Use above function to get ID - gene name conversion table
  tomtom_names <- get_tomtom_names(motif_db_path)
  
  # Read tomtom output table
  all_tom_table <- read.delim(tomtom_path, sep = "\t", quote = "", header = T, 
                              comment.char = "#")
  
  # Convert IDs to gene names
  all_tom_table$Target_ID <- tomtom_names[all_tom_table$Target_ID]
  
  # Remove "(" and ")"
  all_tom_table$Target_ID <- sapply(all_tom_table$Target_ID, function(x){
    str_remove(str_split(x, "\\)")[[1]][1], "\\(")
  })
  
  # Only use alignments below alignment threshold
  all_tom_table <- all_tom_table[all_tom_table$q.value < alignment_threshold,]
  
  return(all_tom_table)
}

scanem_dataset_to_sce <- function(data_path, colData_path){
  require(SingleCellExperiment)
  
  # Read data
  data <- read.csv(data_path, sep="\t")
  
  # Go from comma-separated values to matrix
  lc <- matrix(nrow=nrow(data), ncol=length(str_split(data$ind[1], ",")[[1]]))
  for(i in 1:nrow(lc)){
    lc[i,] <- str_split(data$ind[i], ",")[[1]]
  }
  lc <- t(apply(lc, 1, as.numeric))
  
  # Read metadata
  data_colData <- read.csv(colData_path, sep="\t")
  rownames(lc) <- data$gene
  colnames(lc) <- rownames(data_colData)
  data_colData <- as.data.frame(data_colData)
  data_rowData <- data.frame(row.names=rownames(lc), feature_symbol=rownames(lc))
  
  # Construct SCE object
  sce <- SingleCellExperiment(assays=list(logcounts=as.matrix(lc)), 
                                   colData=as.data.frame(data_colData), 
                                   rowData=as.data.frame(data_rowData))
  # Remove NA values
  sce <- sce[,!is.na(sce$cell_type1)]
}

filter_tomtom_table <- function(tomtom_table, sce){
  expressed_genes <- rownames(sce)[rowSums(logcounts(sce)) > 0]
  tomtom_table <- tomtom_table[tomtom_table$Target_ID %in% expressed_genes,]
  return(tomtom_table)
}

read_network_hdf5 <- function(exp_folder, sce){
  experiment_hdf5_path <- list.files(exp_folder, 
                                     pattern = "M_w\\.h5$", full.names = T)[1]
  
  q <- h5ls(experiment_hdf5_path)
  prefixes <- unique(sapply(q$name, function(x) {str_split(x, "_")[[1]][2] } ))
  
  all_M <- list()
  all_w <- list()
  for(i in 1:length(prefixes)){
    curr_prefix <- prefixes[i]
    
    curr_M <- paste0("M_",curr_prefix)
    curr_M <- rhdf5::h5read(experiment_hdf5_path, curr_M)
    all_M[[curr_prefix]] <- curr_M
    
    curr_w <- paste0("w_",curr_prefix)
    curr_w <- rhdf5::h5read(experiment_hdf5_path, curr_w)
    all_w[[curr_prefix]] <- curr_w
  }
  h5closeAll()
  
  d <- dim(all_M[[1]])[3] # amount of motifs
  m <- dim(all_M[[1]])[2] # length of motif
  
  M_max_activations <- lapply(all_M, get_max_motif_activations)
  
  for(i in 1:length(prefixes)){
    curr_prefix <- prefixes[i]
    
    curr_max_act <- M_max_activations[[curr_prefix]]
    names(curr_max_act) <- paste0(curr_prefix, "_", 0:(d-1))
    
    M_max_activations[[curr_prefix]] <- curr_max_act
  }
  
  scaled_w <- list()
  for(i in 1:length(prefixes)){
    curr_prefix <- prefixes[i]
    
    curr_max_act <- M_max_activations[[curr_prefix]]
    curr_w <- all_w[[curr_prefix]]
    
    # Scale to maximum motif activation
    for(j in 1:nrow(curr_w)){
      curr_w[j,] <- curr_w[j,] * curr_max_act[j]
    }
    
    rownames(curr_w) <- paste0(curr_prefix, "_", 0:(d-1))
    
    scaled_w[[curr_prefix]] <- curr_w
  }
  
  all_w_mat <- do.call(rbind, scaled_w)
  
  colnames(all_w_mat) <- colnames(sce)
  
  return(all_w_mat)
}

create_alignment_graph <- function(tomtom_table){
  require(igraph)
  
  # Get edge connections from table
  edge_vect <- c()
  for(i in 1:nrow(tomtom_table)){
    edge_vect <- c(edge_vect, tomtom_table[i,]$Query_ID, 
                   tomtom_table[i,]$Target_ID)
  }
  edge_weights <- tomtom_table$q.value
  
  # Construct graph
  g <- graph(edge_vect, directed = F)
  
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
  E(g)$arrow.size <- 0.5
  V(g)$size2 <- 3
  
  # Cluster graph
  g <- simplify(g)
  
  return(g)
}

plot_alignment_graph <- function(outdir, alignment_graph, 
                                 alignment_graph_clusters){
  # Plot with clusters
  plot_output_name <- paste0(outdir, "/Motif_alignment_plot_with_clusters.png")
  png(plot_output_name, width = 2000, height = 2000)
  plot(alignment_graph_clusters, alignment_graph, 
       layout = layout_with_graphopt, edge.width=E(alignment_graph)$weight)
  dev.off()
  
  # Plot without clusters
  plot_output_name <- paste0(outdir, 
                             "/Motif_alignment_plot_without_clusters.png")
  png(plot_output_name, width = 2000, height = 2000)
  plot(alignment_graph, 
       layout = layout_with_graphopt, 
       edge.width=E(alignment_graph)$weight)
  dev.off()
}

construct_cluster_df <- function(exp_folder, 
                                 alignment_graph_cluster_groups){
  n_candidates <- sum(str_detect(list.files(exp_folder), "\\.tsv\\.gz"))
  group_reproducibilities <- sapply(alignment_graph_cluster_groups,
                                    FUN=function(x, ncan=n_candidates){
    prefixes <- stringr::str_extract(x, "^best|(^[0-9]+)")
    prefixes <- prefixes[!is.na(prefixes)]
    num_exp_found <- length(unique(prefixes))
    return( (num_exp_found/ncan) )
  })
  cluster_motif_detectors <- sapply(alignment_graph_cluster_groups, 
                                    FUN=function(x, 
                                                 query_options=
                                                   tomtom_table$Query_ID){
    return(x[x %in% query_options])
  })
  cluster_not_motif_detectors <- sapply(alignment_graph_cluster_groups, 
                                        FUN=function(x, 
                                                     query_options=
                                                       tomtom_table$Query_ID){
    return(x[!(x %in% query_options)])
  })
  cluster_motif_detectors <- sapply(X = cluster_motif_detectors, 
                                    FUN=function(X){
    return(paste0(X, collapse=" "))
  })
  cluster_alignments <- sapply(X = cluster_not_motif_detectors, 
                               FUN=function(X){
    return(paste0(X, collapse=" "))
  })
  cluster_df <- data.frame(cluster_motif_detectors, cluster_alignments, 
                           cluster_reproducibility = group_reproducibilities)
  cluster_df <- cluster_df[order(cluster_df$cluster_reproducibility, 
                                 decreasing = T),]
  
  cluster_names <- sapply(cluster_df$cluster_alignments, FUN=function(x){
    if(str_detect(x, "\\(")){
      if(length(str_split(x, " ")[[1]]) > 1){
        if(length(str_split(x, " ")[[1]]) == sum(str_detect(str_split(x, 
                                                                      " ")[[1]], 
                                                            "\\("))){
          curr_names <- str_remove_all(sapply(
            str_split(str_split(x, " ")[[1]], "_"), 
            FUN=function(x) { return(x[1]) }), "\\(|\\)")
          if(length(curr_names) <= 5){
            names <- curr_names
          } else {
            names <- c(curr_names[1:5], "...")
          }
        } else {
          curr_names <- str_split(x, " ")[[1]]
          curr_names <- str_remove_all(sapply(
            str_split(str_split(x, " ")[[1]], "_"), 
            FUN=function(x) { return(x[1]) }), "\\(|\\)")
          if(length(curr_names) <= 5){
            names <- curr_names
          } else {
            names <- c(curr_names[1:5], "...")
          }
        }
      } else {
        curr_names <- str_split(x, " ")[[1]]
        names <- curr_names
      }
    } else {
      curr_names <- str_split(x, " ")[[1]]
      if(length(curr_names) <= 5){
        names <- curr_names
      } else {
        names <- c(curr_names[1:5], "...")
      }
    }
    return(paste(names, collapse = "/"))
  })
  cluster_df$cluster_annot <- cluster_names
  
  return(cluster_df)
}

construct_average_weights_matrix <- function(weights_matrix){
  require(stringr)
  
  celltypes <- unique(str_remove(colnames(weights_matrix), "_pool[0-9]+$"))
  average_weights_matrix <- matrix(nrow=nrow(weights_matrix), 
                                   ncol=length(celltypes))
  for(i in 1:length(celltypes)){
    curr_w_mat <- weights_matrix[,str_detect(colnames(weights_matrix), 
                                             fixed(celltypes[i])), drop=FALSE]
    average_weights_matrix[,i] <- apply(curr_w_mat, 1, mean)
  }
  colnames(average_weights_matrix) <- celltypes
  rownames(average_weights_matrix) <- rownames(weights_matrix)
  average_weights_matrix <- as.data.frame(average_weights_matrix)
  
  return(average_weights_matrix)
}

get_motif_stats <- function(cluster_df,
                            average_weights_matrix, 
                            pseudocount = 1e-10){
  require(mixtools)
  require(ggplot2)
  require(stringr)
  
  # Uses function defined elsewhere in this file
  motif_celltype_entropies <- apply(average_weights_matrix, 1, shan_ent)
  
  num_celltypes <- ncol(average_weights_matrix)
  theoretical_max <- -log2(1/num_celltypes)
  motif_celltype_entropies <- motif_celltype_entropies / theoretical_max
  
  motif_impacts <- apply(average_weights_matrix, 1, function(x) sum(abs(x)))
  
  motif_stats <- data.frame(row.names = names(motif_celltype_entropies), 
                            motif_celltype_entropy=motif_celltype_entropies, 
                            motif_impact=motif_impacts)
  
  mix_model <- normalmixEM(log10(pseudocount+
                                   motif_stats$motif_impact), k=2)
  
  motif_stats$post_1 <- mix_model$posterior[,1]
  motif_stats$post_2 <- mix_model$posterior[,2]
  
  if((sum(motif_stats$post_1 > 0.6) + sum(motif_stats$post_2 > 0.6)) < 
     nrow(motif_stats)){
    cat("No clear bimodal distribution: likely not a lot of 'dead motifs'\n")
    
    good_motifs <- rownames(motif_stats)
  } else {
    cat("Some 'dead motifs' found\n")
    
    good_cluster <- which(mix_model$mu == max(mix_model$mu))
    good_motifs <- rownames(motif_stats)[mix_model$posterior[,good_cluster] > 
                                           0.9]
  }
  
  motif_stats$is_good_motif <- ifelse(rownames(motif_stats) %in% good_motifs, 
                                      "Yes", "No")
  motif_stats$model <- str_remove(rownames(motif_stats), "_[0-9]+$")
  
  q <- ggplot(motif_stats, aes(x = motif_celltype_entropy, 
                          y = motif_impact, 
                          color = is_good_motif)) + 
    geom_point() + scale_x_log10() + scale_y_log10() + theme_bw() + 
    labs(x="Motif 'cell type entropy'", y="Motif impact", color="Good motif?")
  plot(q)
  
  # Get additional stats/annotation
  motif_cluster <- c()
  cluster_df_motifs <- stringr::str_split(cluster_df$cluster_motif_detectors, 
                                          " ")
  for(i in 1:nrow(motif_stats)){
    found <- F
    for(j in 1:length(cluster_df_motifs)){
      if(rownames(motif_stats)[i] %in% cluster_df_motifs[[j]]){
        motif_cluster <- c(motif_cluster, rownames(cluster_df)[j]) 
        found <- T
        break
      }
    }
    if(!found){
      motif_cluster <- c(motif_cluster, "not aligned")
    }
  }
  motif_cluster <- factor(motif_cluster, levels = c("not aligned", 
                                                    rownames(cluster_df)))
  motif_stats$motif_cluster <- motif_cluster
  
  cluster_reprod <- c()
  for(i in 1:nrow(motif_stats)){
    if(motif_stats$motif_cluster[i] == "not aligned") {
      cluster_reprod <- c(cluster_reprod, "NA")
      next
    }
    cluster_reprod <- c(cluster_reprod, 
                        cluster_df[as.character(motif_stats$motif_cluster[i]),
                                   ]$cluster_reproducibility)
  }
  motif_stats$cluster_reprod <- cluster_reprod
  
  cluster_annot <- c()
  for(i in 1:nrow(motif_stats)){
    if(motif_stats$motif_cluster[i] == "not aligned") {
      cluster_annot <- c(cluster_annot, "NA")
      next
    }
    cluster_annot <- c(cluster_annot, 
                       cluster_df[as.character(motif_stats$motif_cluster[i]),
                                  ]$cluster_annot)
  }
  motif_stats$cluster_annot <- cluster_annot
  
  return(motif_stats)
}

get_annotation_colors <- function(motif_stats){
  ann_colors = list(
    cluster_annot = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506",
                      "#a6cee3","#fb9a99","#984ea3","#ffff33", "#040B99")[
                        1:length(unique(motif_stats$cluster_annot))]
  )
  names(ann_colors$cluster_annot) <- names(table(motif_stats$cluster_annot)[
    order(table(motif_stats$cluster_annot), decreasing = TRUE)])
  return(ann_colors)
}

plot_candidate_tf_information <- function(cluster_df, 
                                          weights_matrix, 
                                          motif_stats, 
                                          sce, 
                                          tomtom_table){
  # This will generate plots in the `outdir` that show the expression
  # pattern for the different cluster candidates, including a Spearman 
  # correlation score for the different cluster candidate TFs. 
  
  # In addition to this, it will generate directories per motif cluster that 
  # contain the different plots that led to these correlation scores. 
  
  require(stringr)
  require(magrittr)
  require(pheatmap)
  require(reshape2)
  
  repr_cl_names <- rownames(cluster_df[cluster_df$cluster_reproducibility >= 
                                         0.5,])
  
  weights_across_pools <- list()
  for(i in 1:length(repr_cl_names)){
    curr_cluster <- repr_cl_names[i]
    weights_across_pools[[curr_cluster]] <- colMeans(weights_matrix[
      motif_stats$motif_cluster == curr_cluster,])
  }
  
  all_clusters_corr_dfs <- list()
  all_cluster_corrs <- list()
  all_cluster_mean_exps <- list()
  for(i in 1:length(repr_cl_names)){
    curr_cluster <- repr_cl_names[i]
    
    curr_cluster_tfs <- sapply(str_split(cluster_df[
      curr_cluster,,drop=F]$cluster_alignments, " ")[[1]], FUN=function(x){
      str_remove_all(str_split(x, "_")[[1]][1], "\\(|\\)")
    })
    
    # Keep TFs that are expressed
    curr_cluster_tfs <- curr_cluster_tfs[curr_cluster_tfs %in% 
                                           rownames(sce[rowSums(logcounts(sce)) 
                                                        > 0,])]
    # Get cluster annotation
    curr_cluster_annot <- motif_stats[
      motif_stats$motif_cluster == curr_cluster,,
      drop = FALSE][1,]$cluster_annot
    
    # Candidate TF expression (in pools)
    curr_cluster_df <- logcounts(sce[curr_cluster_tfs,])
    curr_weights_across_pools <- weights_across_pools[[curr_cluster]]
    
    # Get tomtom table for current cluster motifs
    curr_tom_table <- tomtom_table[tomtom_table$Query_ID %in% 
                                     rownames(motif_stats)[
                                       motif_stats$motif_cluster == 
                                         curr_cluster],]
    
    # Sub-directory for cluster correlations
    cluster_outdir <- paste0(outdir, "/Cluster_", 
                             str_remove(
                               str_remove_all(
                                 str_replace_all(
                                   curr_cluster_annot, 
                                   pattern = "/", replacement = "_"), 
                                 "\\."), 
                               "_$"))
    dir.create(cluster_outdir)
    
    if(nrow(curr_cluster_df) != 0){
      
      # Retain TFs that have different expression across pools
      curr_cluster_df <- curr_cluster_df[rowSds(curr_cluster_df) != 0,,
                                         drop=FALSE]
      curr_cluster_df <- curr_cluster_df[,names(curr_weights_across_pools), 
                                         drop=FALSE]
      
      curr_corrs <- c()
      curr_clusters_corr_dfs <- list()
      if(!is.null(nrow(curr_cluster_df))){
        for(j in 1:nrow(curr_cluster_df)){ # Cycle through TFs
          # For the candidate TFs, calculate the correlation between the TF 
          # expression across pools and the weights across pools (Spearman) 
          
          if(sd(curr_cluster_df[j,]) != 0){
            corr <- cor(curr_weights_across_pools, curr_cluster_df[j,], 
                        method="spearman")
            curr_corrs <- c(curr_corrs, corr)
          } else {
            corr <- 0
            curr_corrs <- c(curr_corrs, corr)  
          }
          
          curr_corr_df <- data.frame(x=curr_weights_across_pools, 
                                     y=curr_cluster_df[j,])
          
          curr_corr_df$cell_type <- sce$cell_type1
          curr_corr_df$cell_type_label <- curr_corr_df$cell_type
          curr_corr_df$cell_type_label[
            duplicated(curr_corr_df$cell_type_label)] <- NA
          
          ggplot(curr_corr_df, aes(x,y, color=cell_type)) + geom_point() +
            theme_bw(base_size=14) + xlab("Average motif weight across pools") + 
            ylab(paste0("Expression of ", rownames(curr_cluster_df)[j], 
                        " across pools")) +
            labs(title = paste0(rownames(curr_cluster_df)[j],
                                " (Spearman R = ", round(curr_corrs[j],
                                                         digits=2), ")"), 
                 color="Cell type") +
            theme(legend.position = "none") +
            geom_label_repel(label=curr_corr_df$cell_type_label, alpha=.5) 
          
          ggsave(filename=paste0(cluster_outdir, 
                                 "/Correlation_with_annot_TF_", rownames(curr_cluster_df)[j], 
                                 ".png"), 
                 width = 8, height=8)
          
          ggplot(curr_corr_df, aes(x,y, color=cell_type)) + geom_point() +
            theme_bw(base_size=14) + xlab("Average motif weight across pools") + 
            ylab(paste0("Expression of ", 
                        rownames(curr_cluster_df)[j], " across pools")) +
            labs(title = paste0(rownames(curr_cluster_df)[j],
                                " (Spearman R = ", round(curr_corrs[j],
                                                         digits=2), ")"), 
                 color="Cell type") +
            theme(legend.position = "none")
          
          ggsave(filename=paste0(cluster_outdir, 
                                 "/Correlation_no_annot_TF_", rownames(curr_cluster_df)[j], 
                                 ".png"), 
                 width = 8, height=8)
          
          curr_clusters_corr_dfs[[rownames(curr_cluster_df)[j]]] <- curr_corr_df
        }
      } else {
        if(sd(curr_cluster_df) != 0){
          corr <- cor(curr_weights_across_pools, curr_cluster_df, 
                      method="spearman")
          curr_corrs <- c(curr_corrs, corr)
        } else {
          corr <- 0
          curr_corrs <- c(curr_corrs, corr)  
        }
        
        ggplot(data.frame(x=curr_weights_across_pools, y=curr_cluster_df), aes(x,y)) + geom_point() +
          theme_bw(base_size=14) + xlab("Average motif weight across pools") + ylab(paste0("Expression of ", rownames(curr_cluster_df), " across cell types")) +
          labs(title = paste0("Average motif weight across pools VS expression of ", rownames(curr_cluster_df)),
               subtitle = paste0("Spearman correlation = ", cor(curr_weights_across_pools, curr_cluster_df, method="spearman")))
        
        ggsave(filename=paste0(cluster_outdir, 
                               "/Correlation_no_annot_TF_", rownames(curr_cluster_df)[j], 
                               ".png"), 
               width = 8, height=8)
        
        curr_clusters_corr_dfs[[rownames(curr_cluster_df)]] <- data.frame(x=curr_weights_across_pools, 
                                                                          y=curr_cluster_df)
      }
      names(curr_corrs) <- rownames(curr_cluster_df)
      
      all_cluster_corrs[[curr_cluster]] <- curr_corrs
      all_clusters_corr_dfs[[curr_cluster]] <- curr_clusters_corr_dfs
      all_cluster_mean_exps[[curr_cluster]] <- rowMeans(curr_cluster_df)
      
      curr_row_annot <- data.frame(row.names = rownames(curr_cluster_df), 
                                   cor=curr_corrs)
      curr_row_annot$num_tom_matches <- sapply(rownames(curr_row_annot), 
                                               FUN=function(x){ 
                                                 sum(tomtom_table$Target_ID 
                                                     == x) })
      
      if(!is.null(nrow(curr_cluster_df))){
        curr_cluster_df <- curr_cluster_df[
          order(curr_row_annot$cor, decreasing = T)
          ,,drop=FALSE]
      }
      
      # Calculate expression means in cell types
      curr_cluster_df_avgs <- matrix(nrow=nrow(curr_cluster_df), 
                                     ncol=length(unique(sce$cell_type1)))
      for(z in 1:ncol(curr_cluster_df_avgs)){
        curr_celltype <- unique(sce$cell_type1)[z]
        # drop=FALSE in the case that one cell type contains one cell
        curr_cluster_df_avgs[,z] <- rowMeans(curr_cluster_df[,sce$cell_type1 == 
                                                               curr_celltype, 
                                                             drop=FALSE])
      }
      colnames(curr_cluster_df_avgs) <- unique(sce$cell_type1)
      rownames(curr_cluster_df_avgs) <- rownames(curr_cluster_df)
      
      cluster_name_edited <- str_remove(
        str_remove_all(
          str_replace_all(
            curr_cluster_annot, 
            pattern = "/", replacement = "_"), 
          "\\."), 
        "_$")
      
      img_name <- paste0(outdir, "/Candidate_TF_expression_heatmap_",
                         cluster_name_edited, ".png")
      pheatmap::pheatmap(curr_cluster_df_avgs, 
                         cluster_rows = F, 
                         cluster_cols = F,
                         cellwidth = 15, cellheight = 18, 
                         color=viridis::viridis(100),
                         angle_col = 45, annotation_row = curr_row_annot, 
                         border_color = NA, 
                         width=18, height=15,
                         filename=img_name)
      
    } else {
      cat("Skipping one cluster, no expressed TFs found for cluster")
    }
  }
  
  cluster_correlations_df <- reshape2::melt(all_cluster_mean_exps) %>% 
    magrittr::set_colnames(c("expression", "cluster"))
  cluster_correlations_df$cluster_annot <- sapply(cluster_correlations_df$cluster, 
                                                   function(x){
    motif_stats$cluster_annot[motif_stats$motif_cluster == as.character(x)][1]
  })
  
  all_corrs_df <- reshape2::melt(all_cluster_corrs) %>% 
    magrittr::set_colnames(c("corr", "cluster"))
  cluster_correlations_df$corr <- all_corrs_df$corr
  
  all_cluster_corrs_df <- stack(all_cluster_corrs)
  cluster_correlations_df$TF <- rownames(all_cluster_corrs_df)
  
  return(cluster_correlations_df)
}

plot_motif_family_correlation_overview <- function(cluster_correlations_df, 
                                                   top_n = 5){
  cluster_correlations_df$tf_labels <- NA
  for(i in 1:length(unique(cluster_correlations_df$cluster_annot))){
    curr_cluster <- unique(cluster_correlations_df$cluster)[i]
    curr_corr_tf_df <- cluster_correlations_df[cluster_correlations_df$cluster == curr_cluster,,drop=FALSE]
    curr_corr_tf_df <- curr_corr_tf_df[order(abs(curr_corr_tf_df$corr), decreasing = TRUE),,drop=FALSE]
    annot_tfs <- curr_corr_tf_df$TF[1:(ifelse(nrow(curr_corr_tf_df) < top_n, nrow(curr_corr_tf_df), top_n))]
    cluster_correlations_df$tf_labels[cluster_correlations_df$cluster == curr_cluster &
                                        cluster_correlations_df$TF %in% annot_tfs] <- 
      cluster_correlations_df$TF[cluster_correlations_df$cluster == curr_cluster &
                                    cluster_correlations_df$TF %in% annot_tfs]
  }
  
  q <- ggplot(cluster_correlations_df, aes(x=cluster_annot, y=corr, color=expression)) +
    geom_jitter(width = 0) + 
    theme_bw(base_size=14) + 
    theme(axis.text.x = element_text(angle=45, hjust=1), panel.border = element_rect(colour = NA),
          axis.line = element_line(colour="black")) + 
    labs(x="Motif cluster name", 
         y="Spearman R", color="Mean TF expression across pools") +
    scale_color_viridis_c() +
    geom_label_repel(label=cluster_correlations_df$tf_labels, box.padding = 0.3)
  
  return(q)
}


plot_motif_family_weights_across_celltypes <- function(weights_matrix, 
                                                       sce, 
                                                       motif_stats, 
                                                       motif_family = NA){
  require(reshape2)
  require(ggplot2)
  
  if(is.na(motif_family)){
    stop("motif_family must be supplied")
  }
  
  w_mat_melted <- melt(weights_matrix) %>% 
    magrittr::set_colnames(c("Motif", "Pool", "Weight"))
  w_mat_melted$Celltype <- colData(sce)[w_mat_melted$Pool,"cell_type1"]
  w_mat_melted$`Motif cluster annotation` <- 
    motif_stats[as.character(w_mat_melted$Motif),]$cluster_annot
  
  q <- ggplot(w_mat_melted[w_mat_melted$`Motif cluster annotation` ==
                        motif_family,], 
         aes(x=reorder(Celltype, Weight, mean), y=Weight)) +
    geom_boxplot() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust=1)) + 
    labs(x = "Cell type", y="Average motif weight", 
         title=motif_family)
  
  return(q) 
}
