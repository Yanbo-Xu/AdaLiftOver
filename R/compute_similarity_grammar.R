#' Compute Similarity Grammar
#'
#' This function computes the sequence grammar similarities between query regions and
#' the corresponding list of target regions.
#'
#' @param gr_query GRanges object representing the query regions.
#' @param gr_target_list GRangesList object representing the target regions for each query region.
#' @param hits_query_gr_list GRangesList object of motif hits for the query regions.
#' @param hits_target_gr_list GRangesList object of motif hits for the target regions.
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
    message("Processing query region: ", i, " / ", length(gr_query))
    gr_query_i <- gr_query[i]
    hits_query_gr_i <- hits_query_gr_list[[i]]
    
    # Map motifs in hits_query_gr_i
    if (length(hits_query_gr_i) > 0) {
      hits_query_gr_mapped <- map_motifs(hits_gr = hits_query_gr_i, motif_mapping = motif_mapping,
                                         from_species = 'mouse', to_species = 'human')
    } else {
      hits_query_gr_mapped <- GRanges()
      mcols(hits_query_gr_mapped)$motif <- character(0)
    }
    
    gr_target_regions_i <- gr_target_list[[i]]  # GRanges 对象，包含多个目标区域
    hits_target_gr_list_i <- hits_target_gr_list[[i]]  # GRangesList，与 gr_target_regions_i 对应
    
    similarities_i <- numeric(length(gr_target_regions_i))
    
    for (j in seq_along(gr_target_regions_i)) {
      gr_target_j <- gr_target_regions_i[j]
      hits_target_gr_j <- hits_target_gr_list_i[[j]]
      
      if (length(hits_target_gr_j) == 0) {
        hits_target_gr_j <- GRanges()
        mcols(hits_target_gr_j)$motif <- character(0)
      }
      
      # 计算相似性
      similarity_ij <- compute_similarity_grammar_flat(
        gr_query = gr_query_i,
        gr_target = gr_target_j,
        hits_query_gr = hits_query_gr_mapped,
        hits_target_gr = hits_target_gr_j,
        motif_mapping = motif_mapping,
        all_motifs = all_motifs,
        grammar_size = grammar_size,
        metric = metric
      )
      
      similarities_i[j] <- similarity_ij
    }
    
    # 将相似性添加到目标区域的元数据中
    mcols(gr_target_regions_i)$grammar <- similarities_i
    
    # 保存结果
    result_list[[i]] <- gr_target_regions_i
  }
  
  # 将结果转换为 GRangesList
  result_gr_list <- GRangesList(result_list)
  
  return(result_gr_list)
}

#' Map Motifs Between Species
#'
#' This function maps motifs from the query species to the target species using a provided mapping.
#'
#' @param hits_gr GRanges object containing motif hits with a 'motif' metadata column.
#' @param motif_mapping Data frame containing the mapping between motifs in different species.
#'                      Should have two columns: 'mouse' and 'human'.
#' @param from_species Character, the species of the motifs in hits_gr ('mouse' or 'human').
#' @param to_species Character, the species to map the motifs to ('mouse' or 'human').
#'
#' @return GRanges object with motifs mapped to the target species.
#' @export
map_motifs <- function(hits_gr, motif_mapping, from_species, to_species) {
  stopifnot(class(hits_gr) == 'GRanges')
  
  # 检查 'motif' 列是否存在
  if (!'motif' %in% colnames(mcols(hits_gr))) {
    stop("The hits_gr object must contain a 'motif' metadata column.")
  }
  
  # 确保 from_species 和 to_species 是 'mouse' 或 'human'
  if (!(from_species %in% c('mouse', 'human')) || !(to_species %in% c('mouse', 'human'))) {
    stop("from_species and to_species must be either 'mouse' or 'human'.")
  }
  
  # 创建基序映射：from_species_motif -> to_species_motif
  motif_mapping <- as.data.table(motif_mapping)
  mapping_col_from <- from_species
  mapping_col_to <- to_species
  
  # 创建映射关系
  motif_map <- setNames(motif_mapping[[mapping_col_to]], motif_mapping[[mapping_col_from]])
  
  # 映射基序名称
  hits_gr$motif_mapped <- motif_map[hits_gr$motif]
  
  # 移除无法映射的基序
  hits_gr <- hits_gr[!is.na(hits_gr$motif_mapped)]
  
  # 更新 'motif' 列
  hits_gr$motif <- hits_gr$motif_mapped
  
  # 移除临时列
  hits_gr$motif_mapped <- NULL
  
  return(hits_gr)
}

#' Compute Grammar Matrix
#'
#' This function builds the presence matrix of motifs for the given regions and hits data.
#'
#' @param gr GRanges object representing the regions.
#' @param hits_gr GRanges object of motif hits.
#' @param all_motifs Character vector of all motif names (human motifs).
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
  
  # Filter hits_gr motifs to those in all_motifs
  hits_gr <- hits_gr[hits_gr$motif %in% all_motifs]
  
  # Fill in the matrix
  for (i in seq_along(query_hits)) {
    region_index <- query_hits[i]
    motif_name <- mcols(hits_gr)$motif[subject_hits[i]]
    motif_index <- which(colnames(presence_matrix) == motif_name)
    if (length(motif_index) > 0) {
      presence_matrix[region_index, motif_index] <- TRUE
    }
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
#'                      Should have two columns: 'mouse' and 'human'.
#' @param all_motifs Character vector of all motif names (human motifs).
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
  
  # Map motifs in hits_query_gr (from mouse to human)
  hits_query_gr_mapped <- map_motifs(hits_gr = hits_query_gr, motif_mapping = motif_mapping,
                                     from_species = 'mouse', to_species = 'human')
  
  # Build presence matrices
  query_mat <- compute_grammar_matrix(gr_query, hits_query_gr_mapped, all_motifs, grammar_size)
  target_mat <- compute_grammar_matrix(gr_target, hits_target_gr, all_motifs, grammar_size)
  
  similarity <- compute_similarity_from_matrix_grammar(query_mat, target_mat, metric)
  return(similarity)
}

