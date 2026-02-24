# =============== 聚类空间回投脚本（单样本×RDS输入，直接使用已有聚类） ===============
# v2025-09-04.readRDS
# - 读取已有 Seurat RDS（包含 seurat_clusters 和 UMAP）
# - 与 Visium barcodes_pos.tsv.gz 做【直接精确匹配】
# - 输出：UMAP、空间映射 PDF/PNG、UMAP 投影 RDS
# =====================================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(png)
  library(dplyr)
  library(ggpubr)   # background_image
})

set.seed(1234)

# ===================== 参数区 =====================
# 输入 RDS（包含聚类）
seurat_rds_path <- "path/to/input_file.rds"

# 空间底图与条形码坐标
png_path         <- 'path/to/input_file.png'
barcode_pos_file <- 'path/to/input_file.tsv.gz'

# 输出目录
out_dir        <- 'path/to/output_dir'
analysis_title <- 'Spatial'

custom_cols <-c("0"="#F56867", "1"="#FEB915", "2"="#C798EE","3"="#59BE86","4"="#7495D3","5"="#D1D1D1",
                       "6"="#6D1A9C","7"="#15821E","8"="#3A84E6","9"="#70e014","10"="#787878",
                       "11"="#DB4C6C","12"="#0430e0","13"="#554236","14"="#AF5F3C","15"="#ff7700","16"="#e00417",
                       "17"="#DAB370","18"="#fcfc05","19"="#268785","20"="#09f9f5","21"="#246b93","22"="#cc8e12",
                       "23"="#d561dd","24"="#c93f00","25"="#ddd53e","26"="#4aef7b","27"="#e86502","28"="#9ed84e",
                       "29"="#39ba30","30"="#6ad157","31"="#8249aa","32"="#99db27","33"="#e07233","34"="#ff523f")


# ===================== 读取 RDS =====================
cat("正在读取 Seurat 对象...\n")
object <- readRDS(seurat_rds_path)
cat("读取完成，细胞数:", ncol(object), " 基因数:", nrow(object), "\n")

out_path <- file.path(out_dir, analysis_title)
if(!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

# ===================== UMAP 绘制 =====================
if (!"umap" %in% Reductions(object)) {
  stop("错误: 当前 RDS 不包含 UMAP 投影，请确认保存的对象已经运行过 RunUMAP()。")
}

p_umap <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.8) +
  scale_color_manual(values = custom_cols) +
  ggtitle("UMAP (from existing RDS)")

ggsave(file.path(out_path, "umap_clusters.png"), p_umap, width = 10, height = 8, dpi = 300)
ggsave(file.path(out_path, "umap_clusters.pdf"), p_umap, width = 10, height = 8, dpi = 300)

# 保存 UMAP 投影 RDS（含坐标 + cluster）
umap_embed <- Embeddings(object, "umap")
umap_df <- data.frame(Barcode = rownames(umap_embed), umap_embed, cluster = object$seurat_clusters)
saveRDS(umap_df, file = file.path(out_path, "umap_projection.rds"))

# ===================== 读取底图与 Visium 坐标 =====================
base_img <- readPNG(png_path)

barcode_pos_raw <- read.table(gzfile(barcode_pos_file), header = FALSE, stringsAsFactors = FALSE)
colnames(barcode_pos_raw) <- c("barcode", "pos_w", "pos_h")
barcode_pos <- barcode_pos_raw %>%
  mutate(
    barcode = as.character(barcode),
    pos_w = as.numeric(pos_w),
    pos_h = as.numeric(pos_h)
  )

cal_zoom_rate <- function(width, height) {
  std_width  <- 1000
  std_height <- std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
  if (std_width / std_height > width / height) width / std_width else height / std_height
}
zoom_scale <- cal_zoom_rate(dim(base_img)[2], dim(base_img)[1])
barcode_pos <- barcode_pos %>% mutate(across(c(pos_w, pos_h), ~ .x * zoom_scale))

# ===================== 合并（精确匹配） =====================
meta <- object@meta.data
meta$Barcode <- colnames(object)
meta$cluster <- as.character(meta$seurat_clusters)

cluster_spatial <- meta %>%
  select(Barcode, cluster) %>%
  inner_join(barcode_pos, by = c("Barcode" = "barcode"))

match_rate <- nrow(cluster_spatial) / ncol(object)
cat(sprintf("与坐标表精确匹配成功：%d / %d (%.1f%%)\n",
            nrow(cluster_spatial), ncol(object), 100*match_rate))

# ===================== 空间绘制 =====================
cluster_spatial_plot <- function(df, title = 'Spatial Mapping (all clusters)', point_size = 1) {
  df$cluster <- factor(df$cluster)
  ggplot(df, aes(x = pos_w, y = (dim(base_img)[1] - pos_h), color = cluster)) +
    background_image(base_img) +
    geom_point(shape = 16, alpha = 1, size = point_size) +
    coord_cartesian(xlim = c(0, dim(base_img)[2]), ylim = c(0, dim(base_img)[1]), expand = FALSE) +
    scale_color_manual(values = custom_cols) +
    theme_void() +
    labs(color = "Cluster", title = title) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    )
}

p <- cluster_spatial_plot(cluster_spatial, point_size = 2.4)

ggsave(filename = file.path(out_path, "cluster_mapping.pdf"), plot = p, width = 12, height = 10, dpi = 600)
ggsave(filename = file.path(out_path, "cluster_mapping.png"), plot = p, width = 12, height = 10, dpi = 300)

# ===================== 输出表格 =====================
write.csv(
  data.frame(Barcode = rownames(object@meta.data), seurat_clusters = object$seurat_clusters),
  file = file.path(out_path, "seurat_clusters_all.csv"), row.names = FALSE
)

write.csv(cluster_spatial, file = file.path(out_path, "cluster_spatial_used.csv"), row.names = FALSE)

cat("✅ RDS → 已有聚类 → 全簇空间投影 完成。输出目录:", out_path, "\n")
