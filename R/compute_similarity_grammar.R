#' Compute Similarity Grammar
#'
#' This package provides utilities for motif mapping, grammar computation, and similarity calculations.
#' The main function, `compute_similarity_grammar`, calculates grammar, order, and distance similarities
#' between query and target genomic regions.
#'
#' @import GenomicRanges
#' @import data.table
#' @import Matrix
#' @import transport
#' @export

# LCS function
lcs_length <- function(x, y) {
  m <- length(x)
  n <- length(y)
  L <- matrix(0, nrow = m+1, ncol = n+1)
  for (i in 1:m) {
    for (j in 1:n) {
      if (x[i] == y[j]) {
        L[i+1, j+1] <- L[i, j] + 1
      } else {
        L[i+1, j+1] <- max(L[i, j+1], L[i+1, j])
      }
    }
  }
  return(L[m+1, n+1])
}

# Order similarity
order_similarity <- function(query_seq, target_seq) {
  if (length(query_seq) == 0 || length(target_seq) == 0) {
    return(0)
  }
  lcs_len <- lcs_length(query_seq, target_seq)
  sim <- lcs_len / length(query_seq)
  return(sim)
}

# Extract distance vector
extract_distance_vector <- function(hits_gr) {
  if (length(hits_gr) < 2) return(numeric(0))
  pos <- start(hits_gr)[order(start(hits_gr))]
  diff(pos)
}

# Distance similarity using Wasserstein distance
distance_similarity <- function(query_dist, target_dist) {
  if (length(query_dist) == 0 && length(target_dist) == 0) {
    return(1)
  }
  if (length(query_dist) == 0 || length(target_dist) == 0) {
    return(0)
  }
  w_dist <- wasserstein1d(query_dist, target_dist)
  sim <- 1 / (1 + 0.1 * w_dist)
  return(sim)
}

# Map motifs
map_motifs <- function(hits_gr, motif_mapping, from_species, to_species) {
  stopifnot(class(hits_gr) == 'GRanges')
  if (!'motif' %in% colnames(mcols(hits_gr))) {
    stop("The hits_gr object must contain a 'motif' metadata column.")
  }
  if (!(from_species %in% c('mouse', 'human')) || !(to_species %in% c('mouse', 'human'))) {
    stop("from_species and to_species must be either 'mouse' or 'human'.")
  }
  motif_mapping <- as.data.table(motif_mapping)
  mapping_col_from <- from_species
  mapping_col_to <- to_species
  motif_map <- setNames(motif_mapping[[mapping_col_to]], motif_mapping[[mapping_col_from]])
  hits_gr$motif_mapped <- motif_map[hits_gr$motif]
  hits_gr <- hits_gr[!is.na(hits_gr$motif_mapped)]
  hits_gr$motif <- hits_gr$motif_mapped
  hits_gr$motif_mapped <- NULL
  return(hits_gr)
}

# Compute grammar matrix
compute_grammar_matrix <- function(gr, hits_gr, all_motifs, grammar_size = 500L) {
  stopifnot(class(gr) == 'GRanges', length(gr) > 0, grammar_size >= 150, class(hits_gr) == 'GRanges')
  grammar_size <- as.integer(grammar_size)
  hits_gr <- hits_gr[hits_gr$motif %in% all_motifs]
  offset <- pmax(grammar_size - width(gr), 0) %/% 2
  gr <- gr + offset
  overlaps <- findOverlaps(gr, hits_gr)
  query_hits <- queryHits(overlaps)
  subject_hits <- subjectHits(overlaps)
  presence_matrix <- Matrix(FALSE, nrow = length(gr), ncol = length(all_motifs), sparse = TRUE)
  rownames(presence_matrix) <- as.character(seq_along(gr))
  colnames(presence_matrix) <- all_motifs
  for (i in seq_along(query_hits)) {
    region_index <- query_hits[i]
    motif_name <- mcols(hits_gr)$motif[subject_hits[i]]
    motif_index <- match(motif_name, colnames(presence_matrix))
    if (!is.na(motif_index)) {
      presence_matrix[region_index, motif_index] <- TRUE
    }
  }
  return(presence_matrix)
}

# Compute similarity from matrix grammar
compute_similarity_from_matrix_grammar <- function(query_mat, target_mat, metric = 'cosine') {
  stopifnot(
    nrow(query_mat) == nrow(target_mat) || nrow(query_mat) == 1 || nrow(target_mat) == 1,
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
  } else {
    similarity <- qt_vec / (q_vec + t_vec - qt_vec)
  }
  return(similarity)
}

# Compute similarity grammar flat
compute_similarity_grammar_flat <- function(gr_query, gr_target, hits_query_gr_mapped, hits_target_gr, all_motifs, grammar_size = 500L, metric = 'cosine') {
  stopifnot(
    class(gr_query) == 'GRanges',
    class(gr_target) == 'GRanges',
    length(gr_query) == 1 || length(gr_target) == 1 || length(gr_query) == length(gr_target),
    metric %in% c('cosine', 'jaccard'),
    grammar_size >= 100
  )
  if (length(gr_query) == 0 || length(gr_target) == 0) {
    return(0)
  }
  query_mat <- compute_grammar_matrix(gr_query, hits_query_gr_mapped, all_motifs, grammar_size)
  target_mat <- compute_grammar_matrix(gr_target, hits_target_gr, all_motifs, grammar_size)
  similarity <- compute_similarity_from_matrix_grammar(query_mat, target_mat, metric)
  return(similarity)
}

# Main function to compute similarity grammar
compute_similarity_grammar <- function(gr_query, gr_target_list, hits_query_gr_list, hits_target_gr_list, motif_mapping, all_motifs, grammar_size = 500L, metric = 'cosine', verbose = TRUE) {
  stopifnot(
    class(gr_query) == 'GRanges',
    class(gr_target_list) %in% c('GRangesList', 'CompressedGRangesList'),
    length(gr_query) == length(gr_target_list),
    class(hits_query_gr_list) %in% c('GRangesList', 'CompressedGRangesList'),
    length(gr_query) == length(hits_query_gr_list),
    length(gr_query) == length(hits_target_gr_list),
    metric %in% c('cosine', 'jaccard')
  )
  if (verbose) {
    message('Computing grammar, order, and distance similarity scores (using transport for distance).')
  }
  result_list <- vector("list", length(gr_query))
  for (i in seq_along(gr_query)) {
    if (verbose) {
      message("Processing query region: ", i, " / ", length(gr_query))
    }
    gr_query_i <- gr_query[i]
    hits_query_gr_i <- hits_query_gr_list[[i]]
    if (length(hits_query_gr_i) > 0) {
      hits_query_gr_mapped <- map_motifs(
        hits_gr = hits_query_gr_i, 
        motif_mapping = motif_mapping,
        from_species = 'mouse', 
        to_species = 'human'
      )
    } else {
      hits_query_gr_mapped <- GRanges()
      mcols(hits_query_gr_mapped)$motif <- character(0)
    }
    query_motifs <- as.character(mcols(hits_query_gr_mapped)$motif)
    if (length(hits_query_gr_mapped) > 0) {
      q_ord_idx <- order(start(hits_query_gr_mapped))
      query_motifs <- query_motifs[q_ord_idx]
      query_distances <- extract_distance_vector(hits_query_gr_mapped[q_ord_idx])
    } else {
      query_distances <- numeric(0)
    }
    gr_target_regions_i <- gr_target_list[[i]]
    hits_target_gr_list_i <- hits_target_gr_list[[i]]
    grammar_sims <- numeric(length(gr_target_regions_i))
    order_sims <- numeric(length(gr_target_regions_i))
    dist_sims <- numeric(length(gr_target_regions_i))
    for (j in seq_along(gr_target_regions_i)) {
      gr_target_j <- gr_target_regions_i[j]
      hits_target_gr_j <- hits_target_gr_list_i[[j]]
      if (length(hits_target_gr_j) == 0) {
        hits_target_gr_j <- GRanges()
        mcols(hits_target_gr_j)$motif <- character(0)
      }
      grammar_sims[j] <- compute_similarity_grammar_flat(
        gr_query = gr_query_i,
        gr_target = gr_target_j,
        hits_query_gr_mapped = hits_query_gr_mapped,
        hits_target_gr = hits_target_gr_j,
        all_motifs = all_motifs,
        grammar_size = grammar_size,
        metric = metric
      )
      target_motifs <- as.character(mcols(hits_target_gr_j)$motif)
      if (length(hits_target_gr_j) > 0) {
        t_ord_idx <- order(start(hits_target_gr_j))
        target_motifs <- target_motifs[t_ord_idx]
      }
      order_sims[j] <- order_similarity(query_motifs, target_motifs)
      target_distances <- numeric(0)
      if (length(hits_target_gr_j) > 1) {
        target_distances <- extract_distance_vector(hits_target_gr_j[t_ord_idx])
      }
      dist_sims[j] <- distance_similarity(query_distances, target_distances)
    }
    mcols(gr_target_regions_i)$grammar <- grammar_sims
    mcols(gr_target_regions_i)$order <- order_sims
    mcols(gr_target_regions_i)$distance <- dist_sims
    result_list[[i]] <- gr_target_regions_i
  }
  result_gr_list <- GRangesList(result_list)
  return(result_gr_list)
}