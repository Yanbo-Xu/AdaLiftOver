#' Filter Top Candidate Regions
#'
#' This function filters candidate target regions based on the top percentile
#' of grammar similarity scores.
#'
#' @param gr_candidate A \code{\link[GenomicRanges]{GRangesList}} object with a 'grammar' meta column.
#' @param best_k Integer, the number of top regions to select for each query region. Default is 1.
#' @param top_percentile Numeric, selects regions with grammar scores in the top percentile
#'        (value between 0 and 1). Default is 0.01 (top 1%).
#' @param random Logical, if TRUE, selects \code{best_k} regions randomly; otherwise, selects
#'        the highest-scoring regions. Default is FALSE.
#' @param verbose Logical, whether to print messages. Default is TRUE.
#'
#' @return A \code{\link[GenomicRanges]{GRangesList}} object with filtered regions.
#' @export
gr_candidate_filter <- function(gr_candidate,
                                best_k = 1L,
                                top_percentile = 0.01,
                                random = FALSE,
                                verbose = TRUE) {
  stopifnot(class(gr_candidate) %in% c('CompressedGRangesList', 'GRangesList'),
            best_k >= 1,
            top_percentile > 0 && top_percentile <= 1)
  
  best_k <- as.integer(best_k)
  
  gr_candidate_dt <- as.data.table(gr_candidate)
  gr_candidate_dt$`group_name` <- NULL
  
  stopifnot('grammar' %in% colnames(gr_candidate_dt))
  
  if (verbose) {
    message('Filtering regions based on top percentile grammar similarity.')
  }
  
  # Replace NA with 0 for grammar scores
  gr_candidate_dt[is.na(gr_candidate_dt)] <- 0
  
  threshold <- quantile(gr_candidate_dt$grammar, probs = 1 - top_percentile, na.rm = TRUE)
  
  if (verbose) {
    message(sprintf("Using grammar similarity threshold: %.4f (Top %.0f%%)",
                    threshold, top_percentile * 100))
  }
  
  gr_candidate_dt <- gr_candidate_dt[grammar >= threshold]
  
  if (random) {
    gr_candidate_dt <- gr_candidate_dt[
      , .SD[sample(1:.N, min(.N, best_k))], by = group
    ]
  } else {
    gr_candidate_dt <- gr_candidate_dt[
      , .SD[order(-grammar)][1:min(.N, best_k)], by = group
    ]
  }
  
  gr_candidate_dt$group <- factor(gr_candidate_dt$group, levels = 1:length(gr_candidate))
  
  gr_candidate <- makeGRangesListFromDataFrame(
    gr_candidate_dt,
    split.field = 'group',
    keep.extra.columns = TRUE
  )
  names(gr_candidate) <- NULL
  
  return(gr_candidate)
}