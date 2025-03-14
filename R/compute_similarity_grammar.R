#' Update the motif group counts during the comparison of query and target regions.
#'
#' @param query_vec A logical vector indicating the presence of motif groups in the query region.
#' @param target_vec A logical vector indicating the presence of motif groups in the target region.
#' @param motif_count A vector where each element tracks the count for a specific motif group.
update_motif_count <- function(query_vec, target_vec, motif_count) {
  # Ensure query_vec and target_vec are logical vectors of equal length
  if (length(query_vec) != length(target_vec)) {
    stop("query_vec and target_vec must be of the same length")
  }
  
  # For each motif group, if both query and target contain the motif, increment the count
  for (i in seq_along(query_vec)) {
    if (query_vec[i] && target_vec[i]) {
      motif_count[i] <- motif_count[i] + 1  # Increment count for the motif group at index i
    }
  }
  return(motif_count)
}


###############################################################################
##                   NEW & MODIFIED CODE FOR motif GROUP LOGIC               ##
###############################################################################

#' Parse a single motif group string to identify AND/OR relationships
#'
#' Examples of motif_str:
#'   - "pattern_19+pattern_23" => AND logic
#'   - "pattern_10/pattern_16" => OR  logic
#'   - "pattern_11" (no delimiter) => single motif, treat as OR with one item
parse_motif_group_string <- function(motif_str) {
  motif_str <- trimws(motif_str)
  
  if (grepl("\\+", motif_str) && grepl("/", motif_str)) {
    stop("Mixed '+' and '/' found in the same group; please handle complex logic manually.")
  } else if (grepl("\\+", motif_str)) {
    # AND 逻辑
    motifs <- unlist(strsplit(motif_str, "\\+"))
    motifs <- trimws(motifs)
    return(list(type = "AND", motifs = motifs))
  } else if (grepl("/", motif_str)) {
    # OR 逻辑
    motifs <- unlist(strsplit(motif_str, "/"))
    motifs <- trimws(motifs)
    return(list(type = "OR", motifs = motifs))
  } else {
    # 没有任何分隔符，视为单个 motif
    return(list(type = "OR", motifs = motif_str))
  }
}

#' Parse a data.frame-based motif_mapping to build group-based logic
#'
#' For each row in motif_mapping, we parse from_col and to_col
#' to get from_group, to_group with AND/OR logic.
#' E.g.:
#'   motif_mapping$mouse = "pattern_19+pattern_23"
#'   motif_mapping$human = "pattern_13"
#' => group_list[[i]] = list(
#'      from_group = list(type="AND", motifs=c("pattern_19","pattern_23")),
#'      to_group   = list(type="OR",  motifs=c("pattern_13"))
#'    )
parse_motif_mapping <- function(motif_mapping, from_col = "from_col", to_col = "to_col") {
  stopifnot(from_col %in% names(motif_mapping),
            to_col   %in% names(motif_mapping))
  n <- nrow(motif_mapping)
  group_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    from_str <- motif_mapping[[from_col]][i]
    to_str   <- motif_mapping[[to_col]][i]
    group_list[[i]] <- list(
      from_group = parse_motif_group_string(from_str),
      to_group   = parse_motif_group_string(to_str)
    )
  }
  return(group_list)
}

#' Evaluate group presence for a single region (Query or Target side)
#'
#' @param region_gr   A single GRanges (length=1) representing the region.
#' @param hits_region A GRanges of motif hits (subset) overlapping region_gr
#' @param group_list  The parsed list of group info from parse_motif_mapping()
#' @param side        'from' or 'to', deciding which group_list[[i]]$from_group or $to_group is used
#' @return A logical vector of length = length(group_list), indicating each group is present or not
evaluate_group_presence_single_region <- function(region_gr,
                                                  hits_region,
                                                  group_list,
                                                  side = c("from","to"),
                                                  grammar_size = 500L) {
  side <- match.arg(side)
  
  if (length(hits_region) == 0) {
    # no hits => all groups = FALSE
    return(rep(FALSE, length(group_list)))
  }
  
  # grammar_size offset 处理(可选)
  # region_gr <- region_gr + ...
  
  # 收集 region 中实际出现的 distinct motif
  motif_in_region <- unique(hits_region$motif)
  
  presence_vec <- logical(length(group_list))
  
  for (g in seq_along(group_list)) {
    group_info <- if (side == "from") {
      group_list[[g]]$from_group
    } else {
      group_list[[g]]$to_group
    }
    group_type   <- group_info$type
    group_motifs <- group_info$motifs
    
    if (group_type == "OR") {
      # 只要任何一个 motif_in_region 中有与 group_motifs 匹配的即可
      presence_vec[g] <- any(group_motifs %in% motif_in_region)
    } else {
      # AND 逻辑
      presence_vec[g] <- all(group_motifs %in% motif_in_region)
    }
  }
  
  return(presence_vec)
}

#' Compute similarity for two boolean vectors (same length)
#'   metric = "cosine" or "jaccard"
similarity_one_pair_of_group_vectors <- function(query_vec, target_vec, metric = c("cosine","jaccard")) {
  metric <- match.arg(metric)
  q_num <- as.integer(query_vec)
  t_num <- as.integer(target_vec)
  
  intersect_ <- sum(q_num & t_num)
  q_sum      <- sum(q_num)
  t_sum      <- sum(t_num)
  
  if (metric == "cosine") {
    if (q_sum == 0 || t_sum == 0) {
      return(0)
    } else {
      return(intersect_ / sqrt(q_sum * t_sum))
    }
  } else {
    # jaccard
    union_ <- q_sum + t_sum - intersect_
    if (union_ == 0) {
      return(0)
    } else {
      return(intersect_ / union_)
    }
  }
}

###############################################################################
##               MODIFIED compute_similarity_grammar() + FLAT()              ##
###############################################################################

#' Compute Similarity Grammar (Group-based)
#'
#' Computes group-based sequence grammar similarities between query regions
#' and the corresponding list of target regions. Replaces the old motif-based
#' logic with AND/OR group logic parsed from \code{motif_mapping}.
#'
#' @param gr_query GRanges object representing the query regions.
#' @param gr_target_list GRangesList object representing the target regions for each query region.
#' @param hits_query_gr_list GRangesList object of motif hits for the query regions (in source species).
#' @param hits_target_gr_list GRangesList object of motif hits for the target regions (in target species).
#' @param motif_mapping Data frame containing the mapping between motifs in different species,
#'                      each row describing a group in 'from_species' and 'to_species' columns
#'                      using OR (`/`) or AND (`+`) logic. 
#' @param all_motifs (Optional) If you still want to filter out some motifs, 
#'                   you can supply a vector of recognized motif names. 
#'                   However, for group-based logic, you might not rely on it.
#' @param from_species (Deprecated in group-based logic) 
#' @param to_species   (Deprecated in group-based logic)
#' @param grammar_size Integer, the region extension for each region if needed.
#' @param metric Character, the similarity metric to use ('cosine' or 'jaccard').
#' @param verbose Logical, whether to print messages.
#'
#' @return A GRangesList object with the same structure as gr_target_list, 
#'         with an added metadata column 'grammar' containing the similarity scores.
#' @export
compute_similarity_grammar <- function(gr_query,
                                       gr_target_list,
                                       hits_query_gr_list,
                                       hits_target_gr_list,
                                       motif_mapping,
                                       grammar_size = 500L,
                                       from_col = "mouse",
                                       to_col   = "human",
                                       anno_from_col = "annotation_1",
                                       anno_to_col   = "annotation_2",
                                       metric = 'cosine',
                                       verbose = TRUE) {
  # Step A: parse motif mapping, init group_count
  group_list <- parse_motif_mapping(
    motif_mapping,
    from_col= from_col,
    to_col= to_col,
    anno_from_col= anno_from_col,
    anno_to_col= anno_to_col
  )
  
  n_groups   <- length(group_list)
  motif_count <- integer(n_groups)
  
  colnames_query <- character(n_groups)
  colnames_target <- character(n_groups)
  for (k in seq_len(n_groups)) {
    colnames_query[k]  <- group_list[[k]]$from_annotation
    colnames_target[k] <- group_list[[k]]$to_annotation
  }
  
  
  # Step B: compute Query presence matrix (one row per query region)
  n_query <- length(gr_query)
  if (verbose) message("Building query boolean matrix, number of query regions = ", n_query)
  
  query_boolean_matrix <- matrix(FALSE, nrow=n_query, ncol=n_groups)
  colnames(query_boolean_matrix) <- colnames_query
  
  for (i in seq_len(n_query)) {
    gr_query_i      <- gr_query[i]
    hits_query_gr_i <- hits_query_gr_list[[i]]
    
    query_vec <- evaluate_group_presence_single_region(
      region_gr   = gr_query_i,
      hits_region = hits_query_gr_i,
      group_list  = group_list,
      side        = "from",
      grammar_size = grammar_size
    )
    query_boolean_matrix[i, ] <- query_vec
  }
  
  # Step C: build global target presence matrix
  # first figure out how many total target regions there are
  target_lengths <- sapply(gr_target_list, length)
  n_target_total <- sum(target_lengths)
  
  if (verbose) message("Number of target regions in total = ", n_target_total)
  
  target_boolean_matrix <- matrix(FALSE, nrow=n_target_total, ncol=n_groups)
  colnames(target_boolean_matrix) <- colnames_target
  
  # optional: keep a data frame that records (i, j) => row index
  # e.g. target_info = data.frame( query_index=..., j_index=..., etc. )
  
  current_tindex <- 1
  
  # Step D: main loop - compute similarity & fill in target presence
  # We'll store similarity in the same structure as original
  result_list <- vector("list", n_query)
  
  for (i in seq_len(n_query)) {
    if (verbose) {
      message("Processing query region i = ", i, " / ", n_query)
    }
    n_target_i <- length(gr_target_list[[i]])
    similarities_i <- numeric(n_target_i)
    
    # get the query presence vector from precomputed matrix
    query_vec <- query_boolean_matrix[i, ]
    
    for (j in seq_len(n_target_i)) {
      gr_target_j      <- gr_target_list[[i]][j]
      hits_target_gr_j <- hits_target_gr_list[[i]][[j]]
      
      # compute target presence vector if needed
      target_vec <- evaluate_group_presence_single_region(
        region_gr   = gr_target_j,
        hits_region = hits_target_gr_j,
        group_list  = group_list,
        side        = "to",
        grammar_size = grammar_size
      )
      
      # store in target_boolean_matrix
      target_boolean_matrix[current_tindex, ] <- target_vec
      
      # compute similarity
      sim_ij <- similarity_one_pair_of_group_vectors(query_vec, target_vec, metric)
      similarities_i[j] <- sim_ij
      
      # update motif_count
      motif_count <- update_motif_count(query_vec, target_vec, motif_count)
      
      current_tindex <- current_tindex + 1
    }
    # write back similarity to GRanges
    gr_target_regions_i <- gr_target_list[[i]]
    mcols(gr_target_regions_i)$grammar <- similarities_i
    result_list[[i]] <- gr_target_regions_i
  }
  
  # build final gr_list
  result_gr_list <- GRangesList(result_list)
  
  motif_count_df <- data.frame(
    from_annotation = sapply(group_list, function(g) g$from_annotation),
    to_annotation   = sapply(group_list, function(g) g$to_annotation),
    count           = motif_count
  )
  
  # Step E: return everything
  return(list(
    gr_list               = result_gr_list,
    motif_count           = motif_count_df,
    query_boolean_matrix  = query_boolean_matrix,
    target_boolean_matrix = target_boolean_matrix
  ))
}

parse_motif_mapping <- function(motif_mapping,
                                from_col = "from_col",
                                to_col   = "to_col",
                                anno_from_col = "annotation_1",
                                anno_to_col   = "annotation_2") {
  
  stopifnot(from_col %in% names(motif_mapping),
            to_col   %in% names(motif_mapping),
            anno_from_col %in% names(motif_mapping),
            anno_to_col   %in% names(motif_mapping))
  
  n <- nrow(motif_mapping)
  group_list <- vector("list", n)
  
  for (i in seq_len(n)) {
    from_str <- motif_mapping[[from_col]][i]
    to_str   <- motif_mapping[[to_col]][i]
    
    # 第三、四列：注释信息
    from_anno_str <- motif_mapping[[anno_from_col]][i]
    to_anno_str   <- motif_mapping[[anno_to_col]][i]
    
    group_list[[i]] <- list(
      # 保持原有逻辑
      from_group = parse_motif_group_string(from_str),
      to_group   = parse_motif_group_string(to_str),
      # 新增注释
      from_annotation = from_anno_str,
      to_annotation   = to_anno_str
    )
  }
  return(group_list)
}


#' compute_similarity_grammar_flat (Group-based)
#'
#' This function computes the sequence grammar similarity between
#' one query region and one target region, using group-based logic.
#'
#' @param gr_query A single GRanges (length=1) for the query region
#' @param gr_target A single GRanges (length=1) for the target region
#' @param hits_query_gr GRanges of motif hits for the query region
#' @param hits_target_gr GRanges of motif hits for the target region
#' @param group_list A list from parse_motif_mapping()
#' @param grammar_size integer
#' @param metric 'cosine' or 'jaccard'
#' @return numeric(1) similarity score, or NULL if either region is empty
#' @export
compute_similarity_grammar_flat <- function(gr_query,
                                           gr_target,
                                           hits_query_gr,
                                           hits_target_gr,
                                           group_list,
                                           grammar_size = 500L,
                                           metric = 'cosine',
                                           motif_count) {
  stopifnot(
    class(gr_query) == 'GRanges',
    class(gr_target) == 'GRanges',
    length(gr_query) <= 1,
    length(gr_target) <= 1,
    metric %in% c('cosine', 'jaccard')
  )
  
  if (length(gr_query) == 0 || length(gr_target) == 0) {
    return(NULL)
  }
  
  # Evaluate group presence for query side
  query_vec <- evaluate_group_presence_single_region(
    region_gr   = gr_query,
    hits_region = hits_query_gr,
    group_list  = group_list,
    side        = "from",
    grammar_size = grammar_size
  )
  
  # Evaluate group presence for target side
  target_vec <- evaluate_group_presence_single_region(
    region_gr   = gr_target,
    hits_region = hits_target_gr,
    group_list  = group_list,
    side        = "to",
    grammar_size = grammar_size
  )
  
  # Update the motif count based on query_vec and target_vec
  motif_count <- update_motif_count(query_vec, target_vec, motif_count)
  
  # Compute similarity
  sim <- similarity_one_pair_of_group_vectors(query_vec, target_vec, metric)
  return(list(similarity = sim, motif_count = motif_count))
}
