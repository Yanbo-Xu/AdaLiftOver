#' generate_hits_query_gr_list
#'
#' This function generates a GRangesList object for the query regions by mapping motif hits.
#'
#' @param hits_query Data frame containing motif hits for the query species.
#'                   The data frame should have columns: 'chr', 'start', 'end', and 'motif_name'.
#'                   Can accept finemo hit calling result: hits.tsv.
#' @param gr_query A GRanges object representing the query regions.
#'
#' @return A GRangesList object where each element contains the motif hits for a corresponding query region.
#' @import GenomicRanges
#' @import data.table
#' @export
generate_hits_query_gr_list <- function(hits_query, gr_query) {
  # Validate inputs
  stopifnot(is.data.frame(hits_query),
            all(c("chr", "start", "end", "motif_name") %in% colnames(hits_query)),
            class(gr_query) == "GRanges")
  
  # Convert hits to a GRanges object
  hits_query_gr <- GRanges(
    seqnames = hits_query$chr,
    ranges = IRanges(start = hits_query$start, end = hits_query$end),
    motif = hits_query$motif_name
  )
  
  # Find overlaps between hits and query regions
  overlaps_query <- findOverlaps(hits_query_gr, gr_query)
  hits_in_query <- hits_query_gr[queryHits(overlaps_query)]
  hit_query_index_query <- subjectHits(overlaps_query)
  
  # Split hits into a GRangesList by query regions
  hits_query_gr_list <- split(hits_in_query, factor(hit_query_index_query, levels = seq_along(gr_query)))
  
  return(hits_query_gr_list)
}