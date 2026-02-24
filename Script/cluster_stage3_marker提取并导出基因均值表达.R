########################################################################
## 基于 cluster_stage3 的 marker 基因提取 + 输入基因列表输出表达量
########################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

rm(list = ls())

## 参数
opt <- list(
  outdir = "path/to/output_dir",
  indir  = "path/to/input_file.rds",
  gene_csv = "path/to/input_file.csv"   # <-- 你的基因表 (gene_id, gene_name)
)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

## 1) 读入 Seurat 对象
obj <- readRDS(opt$indir)

if(!"cluster_stage3" %in% colnames(obj@meta.data)){
  stop("❌ metadata 中没有 cluster_stage3 列，请检查 Seurat 对象！")
}
Idents(obj) <- "cluster_stage3"

## 2) 设置 assay
if ("SCT" %in% Assays(obj)) {
  DefaultAssay(obj) <- "SCT"
  obj <- PrepSCTFindMarkers(obj)
} else {
  DefaultAssay(obj) <- "RNA"
}

## 3) 读入基因列表
gene_list <- read_csv(opt$gene_csv)

if(!all(c("gene_id","gene_name") %in% colnames(gene_list))){
  stop("❌ 输入 CSV 必须包含 gene_id 和 gene_name 两列！")
}

## 4) 差异分析（可保留完整结果）
markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

## 5) 计算每个 cluster 的平均表达量
avg_exp <- AverageExpression(obj, assays = DefaultAssay(obj), slot = "data")
avg_exp <- avg_exp[[DefaultAssay(obj)]]
avg_exp_df <- as.data.frame(avg_exp)
avg_exp_df$gene_id <- rownames(avg_exp_df)

## 6) 合并 gene_name
res <- gene_list %>%
  left_join(avg_exp_df, by = "gene_id")

## 7) 保存结果
out_file <- file.path(opt$outdir, "gene_list_with_avgexp.csv")
write_csv(res, out_file)

message("✅ 基因列表 + cluster 表达量导出完成！\n",
        "输出文件: ", out_file)
