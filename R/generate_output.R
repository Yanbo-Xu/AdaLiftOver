generate_output <- function(mapping_result,
                            outdir){
  
  gr_list <- mapping_result$gr_list
  motif_count <- mapping_result$motif_count
  query_boolean_matrix <- mapping_result$query_boolean_matrix
  target_boolean_matrix <- mapping_result$target_boolean_matrix
  
  # 1. gr_list
  combined_gr_list <- unlist(gr_list, use.names = FALSE)
  expanded_names <- rep(mcols(gr)$name, elementNROWS(gr_list))
  mcols(combined_gr_list)$name <- expanded_names
  
  df_gr_list <- as.data.frame(combined_gr_list)
  df_gr_list <- df_gr_list[, c("seqnames", "start", "end", "grammar", "name")]
  colnames(df_gr_list) <- c("chr", "start", "end", "score", "query_region")
  
  # 2. motif_count
  write.table(motif_count, paste0(outdir, "/motif_counts.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # 3. query_boolean_matrix
  df_gr <- as.data.frame(gr)
  df_gr <- df_gr[, c("seqnames", "start", "end", "name")]
  colnames(df_gr) <- c("chr", "start", "end", "query_region")
  
  df_gr$TFBS <- character(nrow(df_gr))
  
  for (i in seq_len(nrow(df_gr))) {
    row_presence <- query_boolean_matrix[i, ]
    
    present_idx <- which(row_presence)
    
    if (length(present_idx) > 0) {
      present_colnames <- colnames(query_boolean_matrix)[present_idx]
      df_gr$TFBS[i] <- paste(present_colnames, collapse = ",")
    } else {
      df_gr$TFBS[i] <- ""
    }
  }
  
  write.table(df_gr, paste0(outdir, "/query_region.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # 4. target_boolean_matrix
  df_gr_list$TFBS <- character(nrow(df_gr_list))
  
  for (i in seq_len(nrow(df_gr_list))) {
    row_presence <- target_boolean_matrix[i, ]
    
    present_idx <- which(row_presence)
    
    if (length(present_idx) > 0) {
      present_colnames <- colnames(target_boolean_matrix)[present_idx]
      df_gr_list$TFBS[i] <- paste(present_colnames, collapse = ",")
    } else {
      df_gr_list$TFBS[i] <- ""
    }
  }
  write.table(df_gr_list, paste0(outdir, "/target_region.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  df_filtered <- df_gr_list[df_gr_list$score != 0, ]
  write.table(df_filtered, paste0(outdir, "/target_region_filterd.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
}
