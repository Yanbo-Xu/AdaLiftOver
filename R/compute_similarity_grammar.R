#' Compute Similarity Grammar
#'
#' This function computes the sequence grammar similarities between query regions and
#' the corresponding list of target regions.
#'
#' @param gr_query GRanges object representing the query regions.
#' @param gr_target_list GRangesList object representing the target regions for each query region.
#' @param hits_query_gr_list GRangesList object of motif hits for the query regions (mouse motifs).
#' @param hits_target_gr_list GRangesList object of motif hits for the target regions (human motifs).
#' @param motif_mapping Data frame containing the mapping between motifs in different species.
#'                      Should have two columns: 'mouse' and 'human'.
#' @param all_motifs Character vector of all motif names (human motifs).
#' @param grammar_size Integer, the size of the region to consider around each region.
#' @param metric Character, the similarity metric to use ('cosine' or 'jaccard').
#' @param verbose Logical, whether to print messages.
#'
#' @return GRangesList object with the same structure as gr_target_list, with an added
#' metadata column 'grammar' containing the similarity scores.
#' @export
compute_similarity_grammar <- function(gr_query,
                                       gr_target_list,
                                       hits_query_gr_list,
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
    class(hits_query_gr_list) %in% c('GRangesList', 'CompressedGRangesList'),
    length(gr_query) == length(hits_query_gr_list),
    length(gr_query) == length(hits_target_gr_list),
    metric %in% c('cosine', 'jaccard')
  )
  
  if (verbose) {
    message('Computing grammar similarity scores.')
  }
  
  result_list <- vector("list", length(gr_query))
  
  for (i in seq_along(gr_query)) {
    if (verbose) {
      message("Processing query region: ", i, " / ", length(gr_query))
    }
    gr_query_i <- gr_query[i]
    hits_query_gr_i <- hits_query_gr_list[[i]]
    
    # Map motifs in hits_query_gr_i from mouse to human
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
    
    gr_target_regions_i <- gr_target_list[[i]]  # GRanges for multiple target regions
    hits_target_gr_list_i <- hits_target_gr_list[[i]]  # GRangesList for target regions hits
    
    similarities_i <- numeric(length(gr_target_regions_i))
    
    for (j in seq_along(gr_target_regions_i)) {
      gr_target_j <- gr_target_regions_i[j]
      hits_target_gr_j <- hits_target_gr_list_i[[j]]
      
      if (length(hits_target_gr_j) == 0) {
        hits_target_gr_j <- GRanges()
        mcols(hits_target_gr_j)$motif <- character(0)
      }
      
      # Compute similarity using already mapped hits_query_gr_mapped
      similarity_ij <- compute_similarity_grammar_flat(
        gr_query = gr_query_i,
        gr_target = gr_target_j,
        hits_query_gr_mapped = hits_query_gr_mapped, # 传入已映射好的查询hits
        hits_target_gr = hits_target_gr_j,
        all_motifs = all_motifs,
        grammar_size = grammar_size,
        metric = metric
      )
      
      similarities_i[j] <- similarity_ij
    }
    
    mcols(gr_target_regions_i)$grammar <- similarities_i
    result_list[[i]] <- gr_target_regions_i
  }
  
  result_gr_list <- GRangesList(result_list)
  return(result_gr_list)
}


#' Map Motifs Between Species
#'
#' @param hits_gr GRanges with 'motif' metadata column.
#' @param motif_mapping Data frame with columns 'mouse' and 'human'.
#' @param from_species 'mouse' or 'human'
#' @param to_species 'mouse' or 'human'
#' @export
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


#' Compute Grammar Matrix
#'
#' @param gr GRanges for regions
#' @param hits_gr GRanges for hits
#' @param all_motifs character vector of human motifs
#' @param grammar_size integer
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
  
  # First filter hits_gr by all_motifs
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


#' Compute Similarity From Matrix Grammar
#'
#' @param query_mat sparse matrix
#' @param target_mat sparse matrix
#' @param metric 'cosine' or 'jaccard'
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
  } else {
    similarity <- qt_vec / (q_vec + t_vec - qt_vec)
  }
  
  return(similarity)
}


#' Compute Similarity Grammar Flat
#'
#' This function computes the sequence grammar similarities between one query region and
#' one target region.
#'
#' @param gr_query GRanges representing a single query region
#' @param gr_target GRanges representing one or multiple target regions combined
#' @param hits_query_gr_mapped GRanges of motif hits for query region, already mapped to human
#' @param hits_target_gr GRanges of motif hits for the target region
#' @param all_motifs character vector of human motifs
#' @param grammar_size integer
#' @param metric 'cosine' or 'jaccard'
#' @export
compute_similarity_grammar_flat <- function(gr_query,
                                            gr_target,
                                            hits_query_gr_mapped,
                                            hits_target_gr,
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
    return(NULL)
  }
  
  query_mat <- compute_grammar_matrix(gr_query, hits_query_gr_mapped, all_motifs, grammar_size)
  target_mat <- compute_grammar_matrix(gr_target, hits_target_gr, all_motifs, grammar_size)
  
  similarity <- compute_similarity_from_matrix_grammar(query_mat, target_mat, metric)
  return(similarity)
}