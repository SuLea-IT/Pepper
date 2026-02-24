###############################################################################
# 指定基因在哪个 cluster：读取基因表 -> 计算按 cluster 的平均表达 -> 给出Top cluster
# 输入：
#   1) Seurat对象 (SCT 或 RNA 均可；建议用合并后的 merged_integrated.rds)
#   2) 基因CSV：path/to/input_file.csv（列名：gene_id,gene_name）
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

## ---------- 路径设置（按需修改） ----------
seurat_rds_path <- "path/to/input_file.rds"
gene_csv_path   <- "path/to/input_file.csv"
output_dir      <- file.path(dirname(seurat_rds_path), "gene_cluster_check")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ---------- 读取 Seurat 对象 ----------
obj <- readRDS(seurat_rds_path)
## ---------- 选择 cluster 列（修正版） ----------
meta_cols <- colnames(obj@meta.data)
cand <- meta_cols[grepl("cluster", meta_cols, ignore.case = TRUE)]

cluster_col <- if ("seurat_clusters" %in% meta_cols) {
  "seurat_clusters"
} else if (length(cand) > 0) {
  cand[1]
} else {
  stop("未找到 cluster 列（如 seurat_clusters）")
}

# 取向量而不是 data.frame
cl_vec <- as.character(obj@meta.data[[cluster_col]])

# 因子顺序：纯数字就按数值升序，否则按字母序
if (all(grepl("^\\d+$", cl_vec))) {
  cl_lev <- sort(unique(cl_vec))
} else {
  cl_lev <- sort(unique(cl_vec), na.last = TRUE)
}
cl_fac <- factor(cl_vec, levels = cl_lev)

# 可选：确保名字与细胞一致（有些版本要求）
names(cl_fac) <- colnames(obj)

# 写回 meta.data（这一行不是必须，但方便后续使用）
obj[[cluster_col]] <- cl_fac

# 关键：给 Idents 传“向量/因子”，不要传 data.frame
Idents(obj) <- cl_fac


## ---------- 选择 assay（修正版） ----------
assay_use <- if ("SCT" %in% names(obj@assays)) {
  "SCT"
} else if ("RNA" %in% names(obj@assays)) {
  "RNA"
} else {
  names(obj@assays)[1]   # 兜底：取第一个 assay
}
DefaultAssay(obj) <- assay_use
message("Assay = ", assay_use)

## ---------- 读取基因列表 ----------
genes_df <- read_csv(gene_csv_path, show_col_types = FALSE)
stopifnot(all(c("gene_id","gene_name") %in% colnames(genes_df)))
target_ids <- unique(genes_df$gene_id)

## ---------- 检查基因是否存在 ----------
all_genes <- rownames(obj[[assay_use]])
found_ids <- intersect(target_ids, all_genes)
missing_ids <- setdiff(target_ids, all_genes)

cat("找到的基因：", length(found_ids), "/", length(target_ids), "\n")
if (length(missing_ids) > 0) {
  cat("未找到（将忽略）：", paste(missing_ids, collapse = ", "), "\n")
}
stopifnot(length(found_ids) > 0)

genes_found_df <- genes_df %>% filter(gene_id %in% found_ids)

## ---------- 计算按 cluster 的平均表达 ----------
# 用 AverageExpression（slot="data"）做群体均值
avg_list <- AverageExpression(obj, assays = assay_use, group.by = cluster_col,
                              slot = "data", verbose = FALSE)
avg_mat  <- avg_list[[assay_use]]  # 行：基因；列：cluster

# 只取关心的基因
avg_sel <- avg_mat[found_ids, , drop = FALSE]

# 输出矩阵（便于你做热图/检查）
avg_out <- cbind(
  gene_id   = rownames(avg_sel),
  gene_name = genes_found_df$gene_name[match(rownames(avg_sel), genes_found_df$gene_id)],
  as.data.frame(avg_sel, check.names = FALSE)
)
write.csv(avg_out, file.path(output_dir, "avg_expression_per_cluster.csv"), row.names = FALSE)

## ---------- 为每个基因找“最主要 cluster” ----------
# 顶峰 cluster + 次峰 + fold-change
which_max   <- apply(avg_sel, 1, which.max)
top_cluster <- colnames(avg_sel)[which_max]
top_value   <- avg_sel[cbind(1:nrow(avg_sel), which_max)]

# 次高
second_value <- apply(avg_sel, 1, function(x){
  if (length(x) == 1) return(NA_real_)
  sort(x, decreasing = TRUE)[min(2, length(x))]
})
second_cluster <- mapply(function(x, top_v){
  if (length(x) == 1) return(NA_character_)
  ord <- order(x, decreasing = TRUE)
  colnames(avg_sel)[ord][min(2, length(ord))]
}, split(avg_sel, row(avg_sel))[[1]], top_value, SIMPLIFY = TRUE, USE.NAMES = FALSE)

fc_top_vs_second <- top_value / pmax(second_value, 1e-9)

# 也给出“显著表达的 cluster 集合”（阈值 = >= 50% 顶峰）
sig_clusters <- apply(avg_sel, 1, function(x){
  keep <- which(x >= 0.5 * max(x, na.rm = TRUE))
  paste(colnames(avg_sel)[keep], collapse = ";")
})

summary_df <- tibble(
  gene_id    = rownames(avg_sel),
  gene_name  = genes_found_df$gene_name[match(rownames(avg_sel), genes_found_df$gene_id)],
  top_cluster      = top_cluster,
  top_expr         = round(top_value, 4),
  second_cluster   = second_cluster,
  second_expr      = round(second_value, 4),
  fc_top_vs_second = round(fc_top_vs_second, 3),
  `clusters_>=50pct_of_top` = sig_clusters
) %>%
  arrange(top_cluster, desc(fc_top_vs_second), desc(top_expr))

write.csv(summary_df,
          file.path(output_dir, "gene_top_cluster_summary.csv"),
          row.names = FALSE)


## ---------- 控制台友好输出 ----------
cat("\n===== Top cluster 摘要（前几行）=====\n")
print(head(summary_df, 15))
cat("\n文件已写入：\n- ", file.path(output_dir, "avg_expression_per_cluster.csv"),
    "\n- ", file.path(output_dir, "gene_top_cluster_summary.csv"), "\n", sep = "")
