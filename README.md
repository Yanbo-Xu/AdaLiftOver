## AdaLiftOver
**AdaLiftOver** is a handy R package for adaptively identifying orthologous regions across different species. For each query region, AdaLiftOver outputs a scored and filtered list of candidate target regions that are most similar to the query region in terms of regulatory information.


## Installation

**AdaLiftOver** can be downloaded and installed in R by: 

```r
## install.packages("devtools")
devtools::install_github("ThomasDCY/AdaLiftOver", build_vignettes = TRUE)
```

If the installation fails, make sure you can install the following R packages:

```r
## data.table
install.packages("data.table")

## Matrix
install.packages("Matrix")

## PRROC
install.packages("PRROC")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## motifmatchr
## See also https://github.com/GreenleafLab/motifmatchr
devtools::install_github("GreenleafLab/motifmatchr")

## rtracklayer
BiocManager::install("rtracklayer")

## GenomicRanges
BiocManager::install("GenomicRanges")

## The BSgenome packages required to run the examples
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

## To build the vignette, we need the following packages
install.packages("rmarkdown")
BiocManager::install("BiocStyle")
```

## A quick start

Download the UCSC chain file from [mm10.hg38.rbest.chain.gz](http://hgdownload.cse.ucsc.edu/goldenpath/mm10/vsHg38/reciprocalBest/mm10.hg38.rbest.chain.gz) and unzip it.

```r
library(AdaLiftOver)

data("data_example")

## load the ENCODE repertoire
data("epigenome_mm10")
data("epigenome_hg38")

## load the UCSC chain file
chain <- rtracklayer::import.chain("mm10.hg38.rbest.chain")

## map the query regions
gr_list <- adaptive_liftover(gr, chain)

## compute epigenome signal similarity
gr_list <- compute_similarity_epigenome(gr, gr_list, epigenome_mm10, epigenome_hg38)

## compute sequence grammar similarity
data("jaspar_pfm_list")
gr_list <- compute_similarity_grammar(gr, gr_list, "mm10", "hg38", jaspar_pfm_list)

## filter target candidate regions
gr_list_filter <- gr_candidate_filter(
    gr_list,
    best_k = 1L,
    threshold = 0.5
)
```

We might need to learn the parameters for a pair of matched epigenome datasets other than the ENCODE repertoire we provide, especially for model organisms other than mice. AdaLiftOver provides a handy training module to estimate the logistic regression parameters and suggest an optimal score threshold without filtering out too many candidate target regions.

```r
data("training_module_example")
training_module(gr_candidate, gr_true)
```



## For a user defined epigenome repertoire

For other model organisms, e.g. rats, we will need to collect orthologous epigenome repertoire from scratch. We hereby illustrate an example.

We can download the UCSC chain file from rat to human [rn6.hg38.rbest.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/vsRn6/reciprocalBest/rn6.hg38.rbest.chain.gz) and then import the chain file with the following code. 

```r
chain <- rtracklayer::import.chain('rn6.hg38.rbest.chain')
```

Make sure to install the **BSgenome** object for rat genome as well.
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn6")
```

These datasets are organized as **GRangesList** objects **epigenome_rn6_test**, **epigenome_hg38_test**
```{r}
data('rat_example')
```

Please refer to the vignette for the code and more details with leave one out cross validation.

After computing the candidate target regions for each rat epigenome peaks, 
we first label the candidate target regions in the human genome as positives if they overlap with the corresponding human epigenome peaks and as negatives otherwise.
Then, we estimate the parameters with the training_module() function.

```{r}
tissues <- names(epigenome_rn6_test)
training_result <- rbindlist(
  lapply(tissues, function(tissue) {
    dt <- training_module(
        gr_candidate_list[[tissue]], 
        epigenome_hg38_test[[tissue]], 
        max_filter_proportion = 0.5, 
        interaction = FALSE
    )
    return(dt)
}))
training_result
```

We can take the averaged logistic regression parameters as default for rat studies.



See the vignette for more information!

```r
browseVignettes("AdaLiftOver")
```


### Reference

**C. Dong**, and **S. Keles**, "AdaLiftOver: High-resolution identification of orthologous regulatory elements with adaptive liftOver".

# Yanbo Xu editing
1. 修改了`compute_similarity_grammar`这个函数。现在可以实现使用finemo/hit calling的结果替换掉原本的方法motifmatchr，生成的结果仍然是：Boolean Matrix，行为region，列为motif。   
2. 添加了`generate_hits_query_gr_list`和`generate_hits_target_gr_list`这两个函数。用以将finemo/hit calling的结果转化成用来比对motif的格式。  
3. 修改了`gr_candidate_filter`这个函数。现在不考虑epigenome similarity，只根据grammar similarity来进行计算。threshold被换为top_percentile，使用前1%作为阈值（先前的阈值被固定为0.5）。
workflow示例操作：
```r
# load query region
NCC_bed <- "/home/xuyanbo/adaliftover/raw_data/Neural_crest.bed"
gr <- import(NCC_bed, format = "BED")

# load the UCSC chain file
chain <- rtracklayer::import.chain("/home/xuyanbo/adaliftover/reference/mm10.hg38.rbest.chain")

# map query regions
gr_list <- adaptive_liftover(gr, chain)

# prepare query hit calling results
hits_query <- fread("/home/xuyanbo/adaliftover/raw_data/mouse_hits_onlypos.tsv")
hits_query_gr_list <- generate_hits_query_gr_list(hits_query, gr)

# prepare target hit calling results
hits_target <- fread("/home/xuyanbo/adaliftover/raw_data/human_hits_onlypos.tsv")
hits_target_gr_list <- generate_hits_taregt_gr_list(hits_target, gr, gr_list)

# compute sequence grammar similarity
motif_mapping <- fread("/home/xuyanbo/adaliftover/output/test/mouse_human_pattern_mapping.tsv", header = TRUE)
all_motifs <- readLines("/home/xuyanbo/adaliftover/raw_data/all_human_motifs.txt")
gr_list <- compute_similarity_grammar(gr, gr_list, hits_query_gr_list, hits_target_gr_list, motif_mapping, all_motifs)

# filter target candidate regions
gr_list_filter <- gr_candidate_filter(
  gr_list,
  best_k = 1L,
  top_percentile = 0.01
)

# export all region
combined_gr <- unlist(gr_target_list_with_similarity, use.names = FALSE)
expanded_names <- rep(mcols(gr)$name, elementNROWS(gr_list))
mcols(combined_gr)$name <- expanded_names
mcols(combined_gr)
export(combined_gr, "/home/xuyanbo/adaliftover/output/test/modify/all_peaks.bed", format = "BED")
df <- as.data.frame(combined_gr)
write.table(df, "/home/xuyanbo/adaliftover/output/test/modify/all_peaks_score.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# export filter region
merged_gr <- unlist(gr_list_filter)
filtered_names <- rep(mcols(gr)$name, elementNROWS(gr_list_filter)) 
mcols(merged_gr)$name <- filtered_names
export(merged_gr, "/home/xuyanbo/adaliftover/output/test/modify_filter_peaks.bed", format = "BED")
df <- as.data.frame(merged_gr)
write.table(df, "/home/xuyanbo/adaliftover/output/test/modify_filter_peaks.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```
