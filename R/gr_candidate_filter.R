#' Filter Candidate Regions with K-means after Selecting Top Region per Query
#'
#' This function first computes a combined score (geometric mean of grammar, order, distance),
#' then for **each query region**, selects the target region with the highest combined score (best_k=1).
#' After that, it applies k-means clustering to these selected regions (one per query),
#' and chooses the cluster with the highest average combined score.
#'
#' @param gr_candidate A \code{\link[GenomicRanges]{GRangesList}} object with 'grammar', 'order', 
#'        and 'distance' meta columns.
#' @param k_clusters Integer, the number of clusters to use for k-means. Default is 2.
#' @param random Logical, if TRUE, for each query region randomly selects one of the top scoring 
#'        regions if there's a tie. Default is FALSE. (Though here we only pick the top 1 by definition)
#' @param verbose Logical, whether to print messages. Default is TRUE.
#'
#' @return A \code{\link[GenomicRanges]{GRangesList}} object that only contains the regions 
#'         from the best-scoring cluster after k-means.
#' @export
gr_candidate_filter <- function(gr_candidate,
                                k_clusters = 2L,
                                random = FALSE,
                                verbose = TRUE) {
  stopifnot(class(gr_candidate) %in% c('CompressedGRangesList', 'GRangesList'),
            k_clusters >= 2)
  
  # Convert to data.table for processing
  gr_candidate_dt <- as.data.table(gr_candidate)
  gr_candidate_dt$`group_name` <- NULL
  
  # Check for required columns
  stopifnot(
    'grammar' %in% colnames(gr_candidate_dt),
    'order' %in% colnames(gr_candidate_dt),
    'distance' %in% colnames(gr_candidate_dt)
  )
  
  if (verbose) {
    message('Selecting top scoring region per query and then applying k-means clustering.')
  }
  
  # Replace NA with 0
  gr_candidate_dt[is.na(gr_candidate_dt)] <- 0
  
  # Calculate combined_score using geometric mean
  gr_candidate_dt[, combined_score := (grammar * order * distance)^(1/3)]
  
  # For each query region (group), select the highest scoring target region
  # best_k = 1, so we just pick the max combined_score per group
  if (random) {
    # If random is TRUE and there's a tie, we can randomly select
    # But since we pick just the top 1, random would only matter if multiple have the exact same score
    # We'll handle ties by ordering and then picking randomly from ties if needed.
    gr_candidate_dt <- gr_candidate_dt[, {
      max_score <- max(combined_score)
      candidates <- .SD[combined_score == max_score]
      if (nrow(candidates) > 1) {
        candidates[sample(1:nrow(candidates), 1)]
      } else {
        candidates
      }
    }, by = group]
  } else {
    # Select the row with the highest combined_score for each group
    gr_candidate_dt <- gr_candidate_dt[, .SD[which.max(combined_score)], by = group]
  }
  
  # Now we have exactly one region per query group
  # Apply k-means on combined_score to cluster them
  if (verbose) {
    message(sprintf("Applying k-means with k=%d on combined_score of selected regions.", k_clusters))
  }
  
  km_res <- kmeans(gr_candidate_dt$combined_score, centers = k_clusters, nstart = 10)
  gr_candidate_dt[, cluster := km_res$cluster]
  
  # Compute cluster averages
  cluster_stats <- gr_candidate_dt[, .(mean_score = mean(combined_score)), by = cluster]
  best_cluster <- cluster_stats[which.max(mean_score), cluster]
  
  if (verbose) {
    message(sprintf("Selected cluster %d with highest mean combined score: %.4f",
                    best_cluster,
                    cluster_stats[cluster == best_cluster]$mean_score))
  }
  
  # Filter by best_cluster
  gr_candidate_dt <- gr_candidate_dt[cluster == best_cluster]
  
  # Retain group structure: 
  # Note: now we only have one region per group that was selected initially, 
  # but after filtering by cluster, some groups may be lost if their region didn't fall into best_cluster.
  # This will reduce the length of resulting GRangesList.
  gr_candidate_dt$group <- factor(gr_candidate_dt$group, levels = unique(gr_candidate_dt$group))
  
  # Convert back to GRangesList
  gr_candidate <- makeGRangesListFromDataFrame(gr_candidate_dt,
                                               split.field = 'group',
                                               keep.extra.columns = TRUE)
  
  names(gr_candidate) <- NULL
  
  return(gr_candidate)
}
