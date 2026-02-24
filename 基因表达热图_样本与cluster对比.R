###############################################################################
# 基因表达热图 – 样本 + Cluster vs 基因（过滤零方差行 & 去除指定样本）
# author: ChatGPT
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
})

## ---------- 环境清理 ----------
rm(list = ls())

## ---------- 1. 路径设置（请按需修改） ----------
seurat_rds_path <- "path/to/input_file.rds"        # 整合后 Seurat 对象
output_dir      <- "path/to/output_dir"           # 结果输出目录
gene_info_csv   <- "path/to/input_file.csv" # gene_id,gene_name

## ---------- 2. 读取 Seurat 对象 ----------
cat("正在读取 Seurat 对象...\n")
object <- readRDS(seurat_rds_path)
cat("Seurat 对象读取成功 ✓\n")

## ---------- 3. 去掉指定样本 ----------
if ("orig.ident" %in% colnames(object@meta.data)) { 
  cat("原始细胞总数:", ncol(object), "\n")
  object <- subset(object, subset = orig.ident != "Fruit_25DPA_CA39_ST")
  cat("过滤后细胞总数:", ncol(object), "\n")
} else {
  stop("meta.data 中未找到 orig.ident 列，无法根据样本名过滤！")
}

## ---------- 4. 基础信息输出 ----------
cat("- 细胞数量:", ncol(object), "\n")
cat("- 基因数量:", nrow(object), "\n")
cat("- 数据层:", names(object), "\n")

## ---------- 5. Cluster 列选择 ----------
cluster_columns <- grep("cluster", colnames(object@meta.data),
                        value = TRUE, ignore.case = TRUE)
if ("harmony_clusters" %in% colnames(object@meta.data)) {
  cluster_col <- "harmony_clusters"
} else if (length(cluster_columns) > 0) {
  cluster_col <- cluster_columns[1]
} else {
  stop("未找到 cluster 信息列")
}
cat("- 使用 cluster 列:", cluster_col, "\n")

## ---------- 6. Assay 选择 ----------
if ("SCT" %in% names(object)) {
  primary_assay <- "SCT"
} else if ("RNA" %in% names(object)) {
  primary_assay <- "RNA"
} else {
  primary_assay <- names(object)[1]   # 兜底：取第一个 assay
}
cat("- 主要使用数据层:", primary_assay, "\n")

## ---------- 7. 创建输出目录 ----------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ---------- 8. 读取基因信息 ----------
gene_info <- read.csv(gene_info_csv, stringsAsFactors = FALSE)
stopifnot(all(c("gene_id", "gene_name") %in% colnames(gene_info)))
target_genes <- gene_info$gene_id

## ---------- 9. 基因存在性检查 ----------
all_genes     <- rownames(object[[primary_assay]])
found_genes   <- intersect(target_genes, all_genes)
missing_genes <- setdiff(target_genes, all_genes)
cat("找到的基因:", length(found_genes), "/", length(target_genes), "\n")
if (length(found_genes) == 0) stop("没有找到任何目标基因！")
found_gene_info <- gene_info[gene_info$gene_id %in% found_genes, ]

## ---------- 10. 组合标签 & 过滤 Cluster ----------
sample_info   <- object@meta.data$orig.ident
cluster_info  <- object@meta.data[[cluster_col]]
sample_cluster <- paste0(sample_info, "_Cluster", cluster_info)

excluded_clusters <- c("1","6", "7", "9", "10")   # 要排除的 cluster
cells_to_exclude  <- cluster_info %in% excluded_clusters
sample_cluster_filtered <- sample_cluster[!cells_to_exclude]
unique_combinations     <- unique(sample_cluster_filtered)

## ---------- 11. 提取表达矩阵 ----------
expr_mat <- GetAssayData(object, assay = primary_assay, slot = "data")
expr_mat <- expr_mat[found_genes, , drop = FALSE]

## ---------- 12. 计算平均表达 ----------
avg_expr <- matrix(0,
                   nrow = length(found_genes),
                   ncol = length(unique_combinations),
                   dimnames = list(found_genes, unique_combinations))
for (i in seq_along(unique_combinations)) {
  comb <- unique_combinations[i]
  idx_in_filtered <- which(sample_cluster_filtered == comb)
  idx_in_original <- which(!cells_to_exclude)[idx_in_filtered]
  avg_expr[, i]   <- if (length(idx_in_original) > 1)
    rowMeans(as.matrix(expr_mat[, idx_in_original]))
  else
    as.numeric(expr_mat[, idx_in_original])
}

## ---------- 13. 注释信息构建 ----------
sample_names  <- sub("_Cluster.*", "", unique_combinations)
cluster_nums  <- sub(".*_Cluster", "", unique_combinations)

cluster2tissue <- c(
  "0" = "Placenta",
  "2" = "Mesocarp1", "3" = "Exocarp",
  "4" = "Mesocarp2", "5" = "Endocarp", "8" = "Gland"
)
tissue_names <- cluster2tissue[cluster_nums]

clean_sample <- gsub("_ST$", "", sub("^Fruit_", "", sample_names))
new_colnames <- make.unique(paste0(clean_sample, "_", tissue_names), sep = "_")
colnames(avg_expr) <- new_colnames

annotation_col <- data.frame(
  Sample = factor(clean_sample),
  Tissue = factor(tissue_names),
  row.names = new_colnames
)
# 颜色
sample_cols <- setNames(brewer.pal(max(3, length(unique(clean_sample))), "Set2"),
                        unique(clean_sample))
tissue_cols <- setNames(brewer.pal(min(length(unique(tissue_names)), 11), "Spectral"),
                        unique(tissue_names))
ann_colors <- list(Sample = sample_cols, Tissue = tissue_cols)


## ---------- 13b. 列排序：先按 DPA，再按样本名，再按指定组织顺序 ----------
order_tissue <- c("Exocarp","Mesocarp1","Mesocarp2","Endocarp","Placenta","Gland")

# 提取 DPA 数字；没有 DPA 的样本记为 NA（会被排在最后）
dpa_num <- suppressWarnings(as.numeric(
  ifelse(grepl("(\\d+)\\s*DPA", clean_sample, ignore.case = TRUE),
         sub(".*?(\\d+)\\s*DPA.*", "\\1", clean_sample, perl = TRUE),
         NA)
))

# 同一 DPA 内的“样本名”键：去掉前缀的 DPA，保留如 CA39/L13 等；
# 没有 DPA 的样本就用其自身名称
sample_key <- ifelse(grepl("^\\d+\\s*DPA", clean_sample, ignore.case = TRUE),
                     sub("^\\d+\\s*DPA_?", "", clean_sample, perl = TRUE),
                     clean_sample)

# 组织层优先级；未知层放到最后且按字母序稳定
tissue_chr  <- as.character(tissue_names)
tissue_rank <- match(tissue_chr, order_tissue)
if (anyNA(tissue_rank)) {
  unknown_idx   <- which(is.na(tissue_rank))
  if (length(unknown_idx) > 0) {
    extras <- sort(unique(tissue_chr[unknown_idx]))
    tissue_rank[unknown_idx] <- length(order_tissue) +
      match(tissue_chr[unknown_idx], extras)
  }
}

# 排序键：DPA ↑ → 样本名 ↑ → 组织层次 ↑ → 原始位置（稳定排序）
col_order <- order(dpa_num, sample_key, tissue_rank, seq_along(dpa_num), na.last = TRUE)

# 应用到矩阵与注释（注意：此时还没有 display_mat，不要排序它）
avg_expr       <- avg_expr[, col_order, drop = FALSE]
annotation_col <- annotation_col[col_order, , drop = FALSE]

## ---------- 14. 行名替换为基因名 ----------
gene2name      <- setNames(found_gene_info$gene_name, found_gene_info$gene_id)
display_mat    <- avg_expr
rownames(display_mat) <- gene2name[rownames(display_mat)]


## ---------- 15. 若为 log 值则反转换 ----------
if (max(display_mat, na.rm = TRUE) < 10) {
  display_mat <- expm1(display_mat)
  cat("检测到 log 值，已转换为原始表达\n")
}

## ---------- 16. 过滤零方差 / NA 行 ----------
row_var   <- apply(display_mat, 1, var, na.rm = TRUE)
bad_genes <- names(row_var)[row_var == 0 | is.na(row_var) | is.infinite(row_var)]
if (length(bad_genes) > 0) {
  cat("剔除零方差或含 NA/Inf 的基因:", length(bad_genes), "个\n")
  display_mat <- display_mat[setdiff(rownames(display_mat), bad_genes), , drop = FALSE]
}

## ---------- 17. 保存原始矩阵 ----------
write.csv(avg_expr, file = file.path(output_dir, "gene_expression_matrix.csv"))

## ---------- 18. 绘制热图 ----------
heat_fun <- function(mat, scale_opt, fname, main_title) {
  p <- pheatmap(
    mat,
    color = colorRampPalette(c("#3288bd", "white", "#d53e4f"))(100),
    border_color = "white",
    cellwidth = 15, cellheight = 15,
    scale = scale_opt,
    cluster_rows = TRUE,
    cluster_cols = FALSE,     # 关闭列聚类，保持自定义顺序   ### NEW
    clustering_method = "ward.D2",
    treeheight_col = 20, treeheight_row = 20,
    legend = TRUE,
    show_rownames = TRUE, show_colnames = TRUE,
    fontsize_row = 10, fontsize_col = 8,
    angle_col = 45,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    main = main_title,
    breaks = seq(-4, 4, length.out = 101)
  )
  pdf(file.path(output_dir, paste0(fname, ".pdf")), width = 26, height = 26)
  print(p); dev.off()
  png(file.path(output_dir, paste0(fname, ".png")), width = 1200, height = 800, res = 150)
  print(p); dev.off()
}

heat_fun(display_mat, "row",
         "gene_expression_heatmap",
         "基因表达热图 (样本+Cluster vs 基因)")

heat_fun(display_mat, "none",
         "gene_expression_heatmap_raw",
         "基因表达热图 (原始值, 样本+Cluster vs 基因)")

## ---------- 19. 统计摘要 ----------
cat("\n表达矩阵统计信息:\n")
cat("- 基因数量(过滤后):", nrow(display_mat), "\n")
cat("- 样本+Cluster 组合数量:", ncol(display_mat), "\n")
cat("- 表达范围:", round(range(display_mat, finite = TRUE), 3), "\n")
cat("- 非零表达比例:",
    round(sum(display_mat > 0) / length(display_mat) * 100, 1), "%\n")
cat("\n输出文件保存在:", output_dir, "\n")
cat("全部完成！\n")
