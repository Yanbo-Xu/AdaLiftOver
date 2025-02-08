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
                                       all_motifs = NULL,   # optional if needed
                                       from_species = NULL, # no longer used
                                       to_species   = NULL, # no longer used
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
    message('Computing grammar similarity scores [group-based logic].')
  }
  
  # Parse motif_mapping => group_list
  group_list <- parse_motif_mapping(motif_mapping, from_col="mouse", to_col="human")
  n_group <- length(group_list)
  if (verbose) {
    message("Parsed motif_mapping into ", n_group, " group(s).")
  }
  
  result_list <- vector("list", length(gr_query))
  
  for (i in seq_along(gr_query)) {
    if (verbose) {
      message("Processing query region: ", i, " / ", length(gr_query))
    }
    gr_query_i      <- gr_query[i]
    hits_query_gr_i <- hits_query_gr_list[[i]]
    
    gr_target_regions_i    <- gr_target_list[[i]]  
    hits_target_gr_list_i  <- hits_target_gr_list[[i]]
    
    similarities_i <- numeric(length(gr_target_regions_i))
    
    for (j in seq_along(gr_target_regions_i)) {
      gr_target_j      <- gr_target_regions_i[j]
      hits_target_gr_j <- hits_target_gr_list_i[[j]]
      
      # Use the new "flat" function that is group-based
      similarity_ij <- compute_similarity_grammar_flat(
        gr_query_i,
        gr_target_j,
        hits_query_gr_i,
        hits_target_gr_j,
        group_list = group_list,
        grammar_size = grammar_size,
        metric = metric
      )
      
      # 防止返回 NULL 的情况
      if (is.null(similarity_ij)) {
        similarity_ij <- 0
      }
      similarities_i[j] <- similarity_ij
    }
    mcols(gr_target_regions_i)$grammar <- similarities_i
    result_list[[i]] <- gr_target_regions_i
  }
  
  result_gr_list <- GRangesList(result_list)
  return(result_gr_list)
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
                                            metric = 'cosine') {
  stopifnot(
    class(gr_query)  == 'GRanges',
    class(gr_target) == 'GRanges',
    length(gr_query) <= 1,
    length(gr_target) <= 1,
    metric %in% c('cosine','jaccard')
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
  
  # Compute similarity
  sim <- similarity_one_pair_of_group_vectors(query_vec, target_vec, metric)
  return(sim)
}