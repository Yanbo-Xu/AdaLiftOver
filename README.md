### Reference

**C. Dong**, and **S. Keles**, "AdaLiftOver: High-resolution identification of orthologous regulatory elements with adaptive liftOver".

# Yanbo Xu editing
1. 修改了`compute_similarity_grammar.R`，输入的`motif list`需要是"pattern_a/pattern_b"或"pattern_a+pattern_b"的形式，分别表示TF family motif和cooperation motif。取消了`all_motif`的输入，将`motif list`的行数（多少个motif group）作为presence matrix的列数。
2. 增加了对`motif group`的计数，在每一组`query region`和`target region`比较的同时计数。

workflow示例操作：
```r
library(data.table)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)

# load query region
NCC_bed <- "/home/xuyanbo/adaliftover/raw_data/Neural_crest.bed"
gr <- import(NCC_bed, format = "BED")
chain <- rtracklayer::import.chain("/home/xuyanbo/adaliftover/reference/mm10.hg38.rbest.chain")

# map query regions
gr_list <- adaptive_liftover(gr, chain)

# prepare query hit calling results
hits_query <- fread("/home/xuyanbo/adaliftover/raw_data/mouse_hits_onlypos.tsv")
hits_query_gr_list <- generate_hits_query_gr_list(hits_query, gr)

# prepare target hit calling results
hits_target <- fread("/home/xuyanbo/adaliftover/raw_data/mouse_to_P4CNCC/hits_p4cncc_onlypos_filtered.tsv")
hits_target_gr_list <- generate_hits_target_gr_list(hits_target, gr, gr_list)

# compute sequence grammar similarity
motif_mapping <- fread("/home/xuyanbo/adaliftover/raw_data/mouse_to_P4CNCC/mouse_p4cncc_motif_list.tsv", header = TRUE)
mapping_result <- compute_similarity_grammar(gr, gr_list, hits_query_gr_list, hits_target_gr_list, motif_mapping)
gr_list <- mapping_result$gr_list
motif_count <- mapping_result$motif_count

gr_list_filter <- gr_candidate_filter(
  gr_list,
  best_k = 1L,
  top_percentile = 0.05
)

combined_gr <- unlist(gr_list, use.names = FALSE)
expanded_names <- rep(mcols(gr)$name, elementNROWS(gr_list))
mcols(combined_gr)$name <- expanded_names
mcols(combined_gr)
# export(combined_gr, "/home/xuyanbo/adaliftover/raw_data/mouse_to_P2CNCC/all_peaks.bed", format = "BED")
df <- as.data.frame(combined_gr)
write.table(df, "/home/xuyanbo/adaliftover/output/mouse_to_P4CNCC/p4cncc_target_regions_score.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

motif_mapping$motif_count <- motif_count
write.table(motif_mapping, "/home/xuyanbo/adaliftover/output/mouse_to_P4CNCC/p4cncc_motif_counts.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

```
