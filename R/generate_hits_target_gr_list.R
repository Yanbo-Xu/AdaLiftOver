#' generate_hits_target_gr_list
#'
#' This function generates a nested list of GRangesList objects for the target regions by mapping motif hits.
#'
#' @param hits_target Data frame containing motif hits for the target species.
#'                    The data frame should have columns: 'chr', 'start', 'end', and 'motif_name'.
#'                    Can accept finemo hit calling result: hits.tsv.
#' @param gr_query A GRanges object representing the query regions.
#' @param gr_target_list A GRangesList object representing the target regions for each query region.
#'
#' @return A list of GRangesList objects where each element corresponds to a query region and contains
#'         GRanges objects for the hits in the associated target regions.
#' @import GenomicRanges
#' @import data.table
#' @export
generate_hits_target_gr_list <- function(hits_target, gr_query, gr_target_list) {
  # Validate inputs
  stopifnot(is.data.frame(hits_target),
            all(c("chr", "start", "end", "motif_name") %in% colnames(hits_target)),
            class(gr_query) == "GRanges",
            class(gr_target_list) %in% c("GRangesList", "CompressedGRangesList"))
  
  # Convert hits to a GRanges object
  hits_target_gr <- GRanges(
    seqnames = hits_target$chr,
    ranges = IRanges(start = hits_target$start, end = hits_target$end),
    motif = hits_target$motif_name
  )
  
  # Initialize list for storing results
  hits_target_gr_list <- vector("list", length(gr_query))
  
  # 创建一个空的 GRanges 对象，用于填充
  empty_gr <- GRanges()
  
  for (i in seq_along(gr_query)) {
    # Get the target regions corresponding to the i-th query region
    target_regions <- gr_target_list[[i]]  # GRanges object with multiple target regions
    
    # If the target regions are empty, store an empty GRangesList and continue
    if (length(target_regions) == 0) {
      # 此时没有目标区域，则 hits_target_gr_list[[i]] 就是空的 GRangesList
      hits_target_gr_list[[i]] <- GRangesList()
      next
    }
    
    # Find overlaps between target regions and hits
    overlaps_target <- findOverlaps(target_regions, hits_target_gr)
    
    # If no overlaps are found, create a GRangesList of empty GRanges for each target region
    if (length(overlaps_target) == 0) {
      # 使用 rep(list(empty_gr), length(target_regions)) 创建一个包含空 GRanges 的列表
      empty_list <- rep(list(empty_gr), length(target_regions))
      hits_target_gr_list[[i]] <- GRangesList(empty_list)
      next
    }
    
    # Extract hits in the target regions
    hits_in_target <- hits_target_gr[subjectHits(overlaps_target)]
    hit_target_indices <- queryHits(overlaps_target)
    
    # Group hits by their corresponding target regions
    hits_per_target <- split(hits_in_target, factor(hit_target_indices, levels = seq_along(target_regions)))
    
    # hits_per_target 是一个 GRangesList，如有必要可以确认一下类型
    # 确保没有 NULL，每个元素应该都是 GRanges，如果不存在命中的目标区域也会是空的 GRanges
    # 一般无需额外处理，因为 split 对 GRanges 会返回 GRangesList，而不存在 NULL。
    
    # Store the grouped hits in the results list
    hits_target_gr_list[[i]] <- hits_per_target
  }
  
  return(hits_target_gr_list)
}