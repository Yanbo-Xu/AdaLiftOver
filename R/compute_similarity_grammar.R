
#' Compute Similarity Grammar
#'
#' This function computes the sequence grammar similarities between query regions and
#' the corresponding list of target regions.
#'
#' @param gr_query GRanges object representing the query regions.
#' @param gr_target_list GRangesList object representing the target regions for each query region.
#' @param hits_query_gr GRanges object of motif hits for the query regions.
#' @param hits_target_gr_list GRangesList object of motif hits for the target regions.
#' @param motif_mapping Data frame containing the mapping between motifs in different species.
#' @param all_motifs Character vector of all motif names.
#' @param grammar_size Integer, the size of the region to consider around each region.
#' @param metric Character, the similarity metric to use ('cosine' or 'jaccard').
#' @param verbose Logical, whether to print messages.
#'
#' @return GRangesList object with the same structure as gr_target_list, with an added
#' metadata column 'grammar' containing the similarity scores.
#' @export
compute_similarity_grammar <- function(gr_query,
                                       gr_target_list,
                                       hits_query_gr,
                                       hits_target_gr_list,
                                       motif_mapping,
                                       all_motifs,
                                       grammar_size = 500L,
                                       metric = 'cosine',
                                       verbose = TRUE) {
  stopifnot(
    class(gr_query) == 'GRanges',
    class(gr_target_list) %in% c('GRangesList', 'CompressedGRangesList'),
    length(gr_query) == length(gr_target_list),
    class(hits_query_gr) == 'GRanges',
    class(hits_target_gr_list) == 'GRangesList',
    length(gr_query) == length(hits_target_gr_list),
    metric %in% c('cosine', 'jaccard')
  )
  
  if (verbose) {
    message('Computing grammar similarity scores.')
  }
  
  # Expand data
  gr_query_flat <- gr_query[rep(seq_along(gr_query), lengths(gr_target_list))]
  gr_target_flat <- unlist(gr_target_list, use.names = FALSE)
  hits_target_gr_flat <- unlist(hits_target_gr_list, use.names = FALSE)
  hits_query_gr_flat <- hits_query_gr[rep(seq_along(hits_query_gr), lengths(gr_target_list))]
  
  # Map motifs in hits_target_gr_flat
  hits_target_gr_flat_mapped <- map_motifs(hits_target_gr_flat, motif_mapping)
  
  # Ensure hits_query_gr_flat motifs are in all_motifs
  hits_query_gr_flat <- hits_query_gr_flat[hits_query_gr_flat$motif %in% all_motifs]
  
  similarity <- compute_similarity_grammar_flat(
    gr_query = gr_query_flat,
    gr_target = gr_target_flat,
    hits_query_gr = hits_query_gr_flat,
    hits_target_gr = hits_target_gr_flat_mapped,
    motif_mapping = motif_mapping,
    all_motifs = all_motifs,
    grammar_size = grammar_size,
    metric = metric
  )
  
  # Combine results
  gr_target_df <- as.data.table(gr_target_list)
  gr_target_df$group <- factor(gr_target_df$group, levels = seq_along(gr_target_list))
  gr_target_df$`group_name` <- NULL
  gr_target_df$grammar <- similarity
  
  gr_target_list <- makeGRangesListFromDataFrame(gr_target_df,
                                                 split.field = 'group',
                                                 keep.extra.columns = TRUE)
  names(gr_target_list) <- NULL
  
  return(gr_target_list)
}

#' Map Motifs Between Species
#'
#' This function maps motifs from one species to another using a provided mapping.
#'
#' @param hits_gr GRanges object containing motif hits with a 'motif' metadata column.
#' @param motif_mapping Data frame containing the mapping between motifs in different species.
#'
#' @return GRanges object with motifs mapped to the target species.
#' @export
map_motifs <- function(hits_gr, motif_mapping) {
  stopifnot(class(hits_gr) == 'GRanges')
  
  # Check if 'motif' column exists
  if (!'motif' %in% colnames(mcols(hits_gr))) {
    stop("The hits_gr object must contain a 'motif' metadata column.")
  }
  
  # Convert mapping to long format
  motif_mapping_long <- melt(
    motif_mapping,
    id.vars = "mouse",
    measure.vars = names(motif_mapping)[-1],
    variable.name = "human_rank",
    value.name = "human_motif"
  )
  motif_mapping_long <- motif_mapping_long[!is.na(human_motif)]
  human_to_mouse_motif <- setNames(motif_mapping_long$mouse, motif_mapping_long$human_motif)
  
  # Map motifs
  hits_gr$motif_mapped <- human_to_mouse_motif[hits_gr$motif]
  
  # Remove unmapped motifs
  hits_gr <- hits_gr[!is.na(hits_gr$motif_mapped)]
  
  # Update 'motif' column
  hits_gr$motif <- hits_gr$motif_mapped
  
  # Remove temporary column
  hits_gr$motif_mapped <- NULL
  
  return(hits_gr)
}

#' Compute Grammar Matrix
#'
#' This function builds the presence matrix of motifs for the given regions and hits data.
#'
#' @param gr GRanges object representing the regions.
#' @param hits_gr GRanges object of motif hits.
#' @param all_motifs Character vector of all motif names.
#' @param grammar_size Integer, the size of the region to consider around each region.
#'
#' @return Sparse logical matrix indicating presence of motifs in regions.
#' @export
compute_grammar_matrix <- function(gr,
                                   hits_gr,
                                   all_motifs,
                                   grammar_size = 500L) {
  stopifnot(class(gr) == 'GRanges',
            length(gr) > 0,
            grammar_size >= 150,
            class(hits_gr) == 'GRanges')
  
  grammar_size <- as.integer(grammar_size)
  
  # Adjust region lengths
  offset <- pmax(grammar_size - width(gr), 0) %/% 2
  gr <- gr + offset
  
  # Find overlaps between regions and hits
  overlaps <- findOverlaps(gr, hits_gr)
  query_hits <- queryHits(overlaps)
  subject_hits <- subjectHits(overlaps)
  
  # Create empty sparse matrix
  presence_matrix <- Matrix(FALSE, nrow = length(gr), ncol = length(all_motifs), sparse = TRUE)
  rownames(presence_matrix) <- as.character(seq_along(gr))
  colnames(presence_matrix) <- all_motifs
  
  # Fill in the matrix
  for (i in seq_along(query_hits)) {
    region_index <- query_hits[i]
    motif_name <- mcols(hits_gr)$motif[subject_hits[i]]
    motif_index <- which(colnames(presence_matrix) == motif_name)
    presence_matrix[region_index, motif_index] <- TRUE
  }
  
  return(presence_matrix)
}

#' Compute Similarity from Grammar Matrices
#'
#' This function computes the similarity between two sets of regions based on their grammar matrices.
#'
#' @param query_mat Sparse matrix of query regions.
#' @param target_mat Sparse matrix of target regions.
#' @param metric Character, the similarity metric to use ('cosine' or 'jaccard').
#'
#' @return Numeric vector of similarity scores.
#' @export
compute_similarity_from_matrix_grammar <- function(query_mat,
                                                   target_mat,
                                                   metric = 'cosine') {
  stopifnot(
    nrow(query_mat) == nrow(target_mat) ||
      nrow(query_mat) == 1 || nrow(target_mat) == 1,
    ncol(query_mat) == ncol(target_mat),
    metric %in% c('cosine', 'jaccard')
  )
  
  if (nrow(query_mat) != nrow(target_mat)) {
    if (nrow(query_mat) == 1) {
      query_mat <- query_mat[rep(1, nrow(target_mat)), , drop = FALSE]
    } else {
      target_mat <- target_mat[rep(1, nrow(query_mat)), , drop = FALSE]
    }
  }
  
  q_vec <- Matrix::rowSums(query_mat)
  t_vec <- Matrix::rowSums(target_mat)
  qt_vec <- Matrix::rowSums(query_mat & target_mat)
  
  if (metric == 'cosine') {
    similarity <- qt_vec / sqrt(q_vec * t_vec)
  } else if (metric == 'jaccard') {
    similarity <- qt_vec / (q_vec + t_vec - qt_vec)
  }
  
  return(similarity)
}

#' Compute Similarity Grammar Flat
#'
#' This function computes the sequence grammar similarities between query regions and
#' the corresponding target regions.
#'
#' @param gr_query GRanges object representing the query regions.
#' @param gr_target GRanges object representing the target regions.
#' @param hits_query_gr GRanges object of motif hits for the query regions.
#' @param hits_target_gr GRanges object of motif hits for the target regions.
#' @param motif_mapping Data frame containing the mapping between motifs in different species.
#' @param all_motifs Character vector of all motif names.
#' @param grammar_size Integer, the size of the region to consider around each region.
#' @param metric Character, the similarity metric to use ('cosine' or 'jaccard').
#'
#' @return Numeric vector of similarity scores.
#' @export
compute_similarity_grammar_flat <- function(gr_query,
                                            gr_target,
                                            hits_query_gr,
                                            hits_target_gr,
                                            motif_mapping,
                                            all_motifs,
                                            grammar_size = 500L,
                                            metric = 'cosine') {
  stopifnot(
    class(gr_query) == 'GRanges',
    class(gr_target) == 'GRanges',
    length(gr_query) == 1 ||
      length(gr_target) == 1 ||
      length(gr_query) == length(gr_target),
    metric %in% c('cosine', 'jaccard'),
    grammar_size >= 100
  )
  
  if (length(gr_query) == 0 || length(gr_target) == 0) {
    return(NULL) # Must return NULL
  }
  
  # Map motifs in hits_target_gr
  hits_target_gr_mapped <- map_motifs(hits_target_gr, motif_mapping)
  
  # Ensure hits_query_gr motifs are in all_motifs
  hits_query_gr <- hits_query_gr[hits_query_gr$motif %in% all_motifs]
  
  # Build presence matrices
  query_mat <- compute_grammar_matrix(gr_query, hits_query_gr, all_motifs, grammar_size)
  target_mat <- compute_grammar_matrix(gr_target, hits_target_gr_mapped, all_motifs, grammar_size)
  
  similarity <- compute_similarity_from_matrix_grammar(query_mat, target_mat, metric)
  return(similarity)
}
