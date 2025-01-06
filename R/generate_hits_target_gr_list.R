#' generate_hits_target_gr_list
#'
#' This function generates a nested list of GRangesList objects for the target regions by mapping motif hits.
#'
#' @param hits_target Data frame containing motif hits for the \strong{target} species.
#'                    Must have columns: 'chr', 'start', 'end', 'motif_name'.
#' @param gr_query A GRanges object representing the query regions.
#' @param gr_target_list A GRangesList object representing the target regions for each query region.
#'
#' @return A list of GRangesList objects where each element corresponds to a query region and contains
#'         GRanges objects for the hits in the associated target regions.
#' @import GenomicRanges
#' @import data.table
#' @export
generate_hits_target_gr_list <- function(hits_target, gr_query, gr_target_list) {
  stopifnot(is.data.frame(hits_target),
            all(c("chr", "start", "end", "motif_name") %in% colnames(hits_target)),
            class(gr_query) == "GRanges",
            class(gr_target_list) %in% c("GRangesList", "CompressedGRangesList"))
  
  hits_target_gr <- GRanges(
    seqnames = hits_target$chr,
    ranges = IRanges(start = hits_target$start, end = hits_target$end),
    motif = hits_target$motif_name
  )
  
  hits_target_gr_list <- vector("list", length(gr_query))
  empty_gr <- GRanges()
  
  for (i in seq_along(gr_query)) {
    target_regions <- gr_target_list[[i]]
    
    if (length(target_regions) == 0) {
      hits_target_gr_list[[i]] <- GRangesList()
      next
    }
    
    overlaps_target <- findOverlaps(target_regions, hits_target_gr)
    
    if (length(overlaps_target) == 0) {
      empty_list <- rep(list(empty_gr), length(target_regions))
      hits_target_gr_list[[i]] <- GRangesList(empty_list)
      next
    }
    
    hits_in_target <- hits_target_gr[subjectHits(overlaps_target)]
    hit_target_indices <- queryHits(overlaps_target)
    
    hits_per_target <- split(
      hits_in_target,
      factor(hit_target_indices, levels = seq_along(target_regions))
    )
    
    hits_target_gr_list[[i]] <- hits_per_target
  }
  
  return(hits_target_gr_list)
}