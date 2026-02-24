# ===================== 基因在哪个 cluster 中表达 =====================
# v2025-09-04

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(readr)
})

set.seed(1234)

# ---------------- 参数 ----------------
rds_final <- "path/to/input_file.rds"
gene_csv  <- "path/to/input_file.csv"
outdir    <- "path/to/output_dir"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---------------- 读取对象和基因列表 ----------------
obj.final <- readRDS(rds_final)
Idents(obj.final) <- obj.final$cluster_stage3

gene_list <- read_csv(gene_csv)
cat("读取到基因数量：", nrow(gene_list), "\n")

# ---------------- 统计表达情况 ----------------
# 获取表达矩阵
expr_mat <- GetAssayData(obj.final, slot = "data")

# 只保留在对象里的基因
genes_in_obj <- intersect(rownames(expr_mat), gene_list$gene_id)
cat("在对象中找到的基因数量：", length(genes_in_obj), "\n")

# 计算每个 cluster 的平均表达
avg_expr <- AverageExpression(obj.final, features = genes_in_obj, assays = DefaultAssay(obj.final))
avg_expr <- avg_expr[[1]] %>% as.data.frame() %>% rownames_to_column("gene_id")

# 合并基因名
avg_expr <- avg_expr %>% left_join(gene_list, by = "gene_id")

# 保存结果
write_csv(avg_expr, file.path(outdir, "gene_cluster_average_expression.csv"))
cat("✅ 保存结果到：", file.path(outdir, "gene_cluster_average_expression.csv"), "\n")

# ---------------- 可视化 ----------------
# DotPlot
p_dot <- DotPlot(obj.final, features = genes_in_obj, group.by = "cluster_stage3") +
  RotatedAxis() + ggtitle("Gene expression across clusters")
ggsave(file.path(outdir, "gene_cluster_dotplot.png"), p_dot, width = 10, height = 6, dpi = 300)

# 单基因 violin plot（逐个保存）
for (g in genes_in_obj) {
  p_vln <- VlnPlot(obj.final, features = g, group.by = "cluster_stage3", pt.size = 0.1) +
    ggtitle(paste0(g, " (", gene_list$gene_name[gene_list$gene_id == g], ")"))
  ggsave(file.path(outdir, paste0("gene_", g, "_violin.png")),
         p_vln, width = 6, height = 4, dpi = 300)
}

cat("✅ 完成基因表达检查；结果目录：", outdir, "\n")
