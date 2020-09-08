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
          curr_names <- str_remove_all(sapply(str_split(str_split(cluster_df$cluster_alignments[6], " ")[[1]], "_"), FUN=function(x) { return(x[1]) }), "\\(|\\)")
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
                                               rowSums(curr_counts[,indx[indx <= ncol(curr_counts)], drop=F])})
        
        colnames(curr_fragmented_counts) <- paste0(paste0(curr_cell_type, "_pool"), 
                                                   1:ncol(curr_fragmented_counts))
        
        rownames(curr_fragmented_barcodes) <- paste0(paste0(curr_cell_type, "_pool"), 
                                                     1:ncol(curr_fragmented_counts))
        
        total_new_counts <- cbind(total_new_counts, curr_fragmented_counts)
        total_new_barcodes <- rbind(total_new_barcodes, curr_fragmented_barcodes)
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
  celltypes_combined <- stringr::str_remove(colnames(total_new_counts), pattern = "_pool[0-9]+$")
  
  combined_sce <- SingleCellExperiment(assays=list(counts=total_new_counts), 
                                       colData=data.frame(cell_type1=celltypes_combined),
                                       rowData=rowData(sce))
  
  combined_sce <- scater::logNormCounts(combined_sce)
  
  dropout_ratio <- get_sce_sparsity(combined_sce)
  cat("\nDropout ratio:", round(dropout_ratio, 3), "\n")
  cat("Cell number:", ncol(combined_sce), "\n")
  cat("Done\n\n")
  
  return(combined_sce)
}
