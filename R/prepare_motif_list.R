
#' prepare_motif_list
#'
#' Prepare a motif list from PFM.meme and hits.tsv files
#'
#' @param meme_file Path to the PFM.meme file containing motif patterns
#' @param hits_file Path to the hits.tsv file containing hits calling results
#' @param motif_length Length of each motif (default: 30)
#' @return A list containing:
#'         - motif_list: A PWMatrixList object for use in Adaliftover
#'         - hits_data: A GRanges object with motif mapping from hits.tsv
#' @import TFBSTools
#' @import GenomicRanges
#' @import data.table
#' @export
prepare_motif_list <- function(meme_file, hits_file, motif_length = 30) {
  # 1. Load PFM.meme file and convert to PWMatrixList
  message("Processing PFM.meme file to PFMatrixList...")
  motif_list <- parse_meme_to_PFMatrixList(meme_file, motif_length)

  # 2. Load hits.tsv and convert to GRanges
  message("Processing hits.tsv file to GRanges...")
  gr_hits <- parse_hits_to_GRanges(hits_file)
  # 3. Return motif_list and hits_data
  return(list(motif_list = motif_list, hits_data = gr_hits))
}


# convert PFM.meme to PWMatrixList
parse_meme_to_PFMatrixList <- function(meme_file, motif_length = 30) {
  # 读取 MEME 文件
  lines <- readLines(meme_file)
  
  motifs <- list()  # 保存所有 motif 的信息
  current_motif <- NULL
  matrix_data <- NULL
  collecting_matrix <- FALSE  # 是否在收集矩阵数据
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # 跳过文件的全局元数据部分
    if (i <= 6) next
    
    # 检测 MOTIF 开头的行
    if (grepl("^MOTIF", line)) {
      # 如果已经在解析一个 motif，保存之前的
      if (!is.null(current_motif) && !is.null(matrix_data)) {
        profile_matrix <- do.call(rbind, matrix_data)
        if (ncol(profile_matrix) != 4) {
          stop("Invalid motif matrix: must have exactly 4 columns (A, C, G, T).")
        }
        motifs[[current_motif]] <- profile_matrix
        matrix_data <- NULL  # 重置矩阵数据
      }
      
      # 开始新 motif 的解析
      current_motif <- strsplit(line, "\\s+")[[1]][2]
      collecting_matrix <- FALSE  # 暂不收集矩阵
      next
    }
    
    # 检测 letter-probability matrix 行
    if (grepl("^letter-probability matrix", line)) {
      collecting_matrix <- TRUE  # 开始收集矩阵
      row_count <- 0  # 初始化计数器
      next
    }
    
    # 如果在收集矩阵数据，并且当前行是矩阵
    if (collecting_matrix && grepl("^\\s*\\d", line)) {
      row_count <- row_count + 1
      values <- as.numeric(strsplit(line, "\\s+")[[1]])
      if (length(values) == 4) {
        if (is.null(matrix_data)) {
          matrix_data <- list(values)
        } else {
          matrix_data <- c(matrix_data, list(values))
        }
      }
      
      # 如果已经读取了 motif_length 行，完成当前 motif
      if (row_count == motif_length) {
        collecting_matrix <- FALSE  # 停止收集矩阵
      }
      next
    }
  }
  
  # 保存最后一个 motif
  if (!is.null(current_motif) && !is.null(matrix_data)) {
    profile_matrix <- do.call(rbind, matrix_data)
    if (ncol(profile_matrix) != 4) {
      stop("Invalid motif matrix: must have exactly 4 columns (A, C, G, T).")
    }
    motifs[[current_motif]] <- profile_matrix
  }
  
  # 转换为 PFMatrixList
  pfm_list <- lapply(names(motifs), function(name) {
    profile_matrix <- motifs[[name]]
    
    # 转置矩阵并添加行名
    profile_matrix <- t(profile_matrix)
    rownames(profile_matrix) <- c("A", "C", "G", "T")
    
    # 转换为计数矩阵
    count_matrix <- round(profile_matrix * 1000)  # 每个值乘以 1000 转为整数
    
    PFMatrix(ID = name,
             name = name,
             matrixClass = "PFM",
             bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
             profileMatrix = count_matrix)  # 使用计数矩阵
    
  })
  
  pfm_list <- do.call(PFMatrixList, pfm_list)
  
  return(pfm_list)
}

# convert hits calling .tsv to GRanges
parse_hits_to_GRanges <- function(hits_file) {
  # 读取 hits.tsv 文件
  hits_data <- fread(hits_file)
  
  # 验证列名
  stopifnot(all(c("chr", "start", "end", "motif_name") %in% colnames(hits_data)))
  
  # 保留 motif_name 包含 "pos_patterns" 的行
  hits_data <- hits_data[grep("^pos_patterns", hits_data$motif_name)]
  
  # 修改 motif_name 列的值
  hits_data[, motif_name := gsub("^pos_patterns\\.", "", motif_name)]
  
  # 转换为 GRanges
  gr_hits <- GRanges(
    seqnames = hits_data$chr,
    ranges = IRanges(start = hits_data$start, end = hits_data$end),
    pattern = hits_data$motif_name
  )
  
  return(gr_hits)
}