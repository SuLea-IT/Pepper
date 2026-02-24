# =============== 聚类空间回投脚本（单样本×矩阵输入，纯重新聚类，绘制所有簇） ===============
# v2025-09-01.solo
# - 直接从 10x 单样本矩阵读取（filtered_feature_bc_matrix）
# - 本脚本内完成 Normalize/SCT → PCA → 邻居图 → 聚类 → UMAP
# - 与 Visium barcodes_pos.tsv.gz 做【直接精确匹配】（不做后缀裁剪/统一）
# - 绘制【所有】簇的空间分布（不再有 retain_clusters / cluster_to_tissue）
# - 输出：UMAP、空间映射 PDF、重新聚类 RDS、簇表 CSV
# =====================================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(png)
  library(dplyr)
  library(RColorBrewer)
  library(ggpubr)   # background_image
})

set.seed(1234)

# ===================== 参数区 =====================
# 单样本矩阵目录（10x 的 filtered_feature_bc_matrix）
matrix_dir <- "path/to/input_dir"

# 空间底图与条形码坐标（同一张切片的 PNG 与 barcodes_pos.tsv.gz）
png_path         <- 'path/to/input_file.png'
barcode_pos_file <- 'path/to/input_file.tsv.gz'

# 重新聚类相关参数
method <- "NFS"   # 可选："SCT" 或 "NFS"
res    <- 0.5      # 分辨率
npca   <- 30       # PCA 维度数

# 输出目录
out_dir        <- 'path/to/output_dir'
analysis_title <- 'cluster_projection_single'


# ===================== 工具函数（仅缩放/配色） =====================
cal_zoom_rate <- function(width, height) {
  std_width  <- 1000
  std_height <- std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
  if (std_width / std_height > width / height) width / std_width else height / std_height
}

make_palette <- function(n) {
  base <- c("#F56867", "#FEB915", brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
  if (n <= length(base)) base[seq_len(n)] else rep(base, length.out = n)
}

# ===================== 读入矩阵 & 重新聚类 =====================
mat   <- Read10X(matrix_dir)
object <- CreateSeuratObject(counts = mat, assay = "RNA")
cat("读取矩阵完成，细胞数:", ncol(object), "基因数:", nrow(object), "
")

if (method == "SCT") {
  cat("使用 SCTransform 归一化...
")
  object <- SCTransform(object, assay = "RNA", verbose = FALSE)
  assay_for_dimred <- "SCT"
} else {
  cat("使用 NormalizeData/FindVariableFeatures/ScaleData...
")
  object <- NormalizeData(object, assay = "RNA")
  object <- FindVariableFeatures(object, assay = "RNA", selection.method = "vst", nfeatures = 2000)
  object <- ScaleData(object, assay = "RNA", features = rownames(object))
  assay_for_dimred <- "RNA"
}

cat("PCA/Neighbors/Clusters/UMAP ...
")
object <- RunPCA(object, assay = assay_for_dimred, npcs = max(50, npca), verbose = FALSE)
object <- FindNeighbors(object, reduction = "pca", dims = 1:npca)
object <- FindClusters(object, resolution = res)
object <- RunUMAP(object, reduction = "pca", dims = 1:npca)

# UMAP 预览
p_umap <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters", pt.size = 1.2) +
  ggtitle(sprintf("UMAP (res=%.2f, dims=%d)", res, npca))

out_path <- file.path(out_dir, analysis_title)
if(!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

ggsave(file.path(out_path, "umap_clusters.png"), p_umap, width = 10, height = 8, dpi = 300)

# ===================== 读取底图与 Visium 坐标（精确匹配） =====================
base_img <- readPNG(png_path)

barcode_pos_raw <- read.table(gzfile(barcode_pos_file), header = FALSE, stringsAsFactors = FALSE)
colnames(barcode_pos_raw) <- c("barcode", "pos_w", "pos_h")
barcode_pos <- barcode_pos_raw %>%
  mutate(
    barcode = as.character(barcode),
    pos_w = as.numeric(pos_w),
    pos_h = as.numeric(pos_h)
  )

zoom_scale <- cal_zoom_rate(dim(base_img)[2], dim(base_img)[1])
barcode_pos <- barcode_pos %>% mutate(across(c(pos_w, pos_h), ~ .x * zoom_scale))

# ===================== 合并（不做条形码后缀裁剪；直接精确匹配） =====================
meta <- object@meta.data
meta$Barcode <- colnames(object)
meta$cluster <- as.character(meta$seurat_clusters)

cluster_spatial <- meta %>%
  select(Barcode, cluster) %>%
  inner_join(barcode_pos, by = c("Barcode" = "barcode"))

match_rate <- nrow(cluster_spatial) / ncol(object)
cat(sprintf("与坐标表精确匹配成功：%d / %d (%.1f%%)
", nrow(cluster_spatial), ncol(object), 100*match_rate))

# 绘图点大小
point_size <- 2
# ===================== 绘图（绘制所有簇） =====================
cluster_spatial_plot <- function(df, title = 'Spatial Mapping (all clusters)', point_size = 5) {
  df$cluster <- factor(df$cluster)
  cols <- setNames(make_palette(length(levels(df$cluster))), levels(df$cluster))
  
  ggplot(df, aes(x = pos_w, y = (dim(base_img)[1] - pos_h), color = cluster)) +
    background_image(base_img) +
    geom_point(shape = 16, alpha = 0.8, size = point_size) +
    coord_cartesian(xlim = c(0, dim(base_img)[2]), ylim = c(0, dim(base_img)[1]), expand = FALSE) +
    scale_color_manual(values = cols) +
    theme_void() +
    labs(color = "Cluster", title = title) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
}

p <- cluster_spatial_plot(cluster_spatial, point_size = point_size)

ggsave(filename = file.path(out_path, "cluster_mapping.pdf"), plot = p, width = 12, height = 10, dpi = 600)

# ===================== 输出对象与表格 =====================
saveRDS(object, file = file.path(out_path, "object_reclustered.rds"))
write.csv(
  data.frame(Barcode = rownames(object@meta.data), seurat_clusters = object$seurat_clusters),
  file = file.path(out_path, "seurat_clusters_all.csv"), row.names = FALSE
)
# markers & 平均表达量表
write.csv(markers, file = file.path(out_path, "marker_gene_onlypos.csv"), row.names = FALSE)
write.csv(avg_expr, file = file.path(out_path, "average_expression.csv"))


# 空间用于绘图的数据（便于复现）
write.csv(cluster_spatial, file = file.path(out_path, "cluster_spatial_used.csv"), row.names = FALSE)
cat("✅ 单样本矩阵 → 重新聚类 → 全簇空间投影 完成。输出目录:", out_path, "
")
