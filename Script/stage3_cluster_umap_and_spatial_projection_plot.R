# ===================== Plot Stage-3 Clusters =====================
# v2025-09-06
# - 输入：object_stage3_relabel.rds
# - 绘制：seurat_clusters 的 UMAP & 空间投影
# ================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(png)
  library(dplyr)
  library(ggpubr)
})

set.seed(1234)

# ---------------- 参数 ----------------
rds_path   <- "path/to/input_file.rds"
png_path   <- "path/to/input_file.png"
#png_path   <- "path/to/input_file.png"
pos_file   <- "path/to/input_file.tsv.gz"
out_dir    <- "path/to/output_dir"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
# ---------------- 调色板 ----------------f67f7e
make_palette <- function(n) {
  base <- c("#c53e7a", "#c798ee", "#85E045", "#22a5da", "#F55BA2", "#72c898",
            "#833cab", "#5796ea", "#e8e24a", "#15821e",
            "#8c8c8c", "#f67f7e","#89a5da")
  if (n <= length(base)) base[seq_len(n)] else rep(base, length.out = n)
}
# base <- c("#f67f7e", "#fec338", "#C798EE", "#59BE86", "#7495D3", "#DB4C6C",
#           "#15821E", "#3A84E6", "#70E014", "#554236", "#AF5F3C", "#FF7700",
#           "#268785", "#20B2AA", "#931635", "#4682B4", "#FFD700", "#87CEFA")
# if (n <= length(base)) base[seq_len(n)] else rep(base, length.out = n)
# }
# ---------------- 新版缩放函数 ----------------
cal_zoom_rate <- function(barcode_pos, img_w, img_h, tol=0.05) {
  max_w <- max(barcode_pos$pos_w, na.rm=TRUE)
  max_h <- max(barcode_pos$pos_h, na.rm=TRUE)
  if (abs(max_w - img_w)/img_w > tol || abs(max_h - img_h)/img_h > tol) {
    zoom_scale <- min(img_w / max_w, img_h / max_h)
    message("⚠️ 坐标已缩放，比例=", round(zoom_scale, 3))
    barcode_pos <- barcode_pos %>%
      mutate(pos_w = pos_w * zoom_scale,
             pos_h = pos_h * zoom_scale)
  } else {
    message("✅ 坐标范围匹配，无需缩放")
  }
  return(barcode_pos)
}
object <- readRDS(rds_path)


if (!"seurat_clusters" %in% colnames(object@meta.data)) {
  stop("❌ RDS 中未找到 seurat_clusters 列，请检查输入对象")
}

# ---------------- 绘制 UMAP ----------------
n_clusters <- length(unique(object$seurat_clusters))
cols <- make_palette(n_clusters)

p_umap <- DimPlot(object, reduction="umap", group.by="seurat_clusters", pt.size=1.8) +
  scale_color_manual(values=cols) +
  ggtitle("UMAP (seurat_clusters)")

ggsave(file.path(out_dir,"umap_seurat_clusters.png"), p_umap, width=10, height=8, dpi=300)
ggsave(file.path(out_dir,"umap_seurat_clusters.pdf"), p_umap, width=10, height=8, dpi=600)

# ---------------- 空间投影 ----------------
if (file.exists(png_path) && file.exists(pos_file)) {
  base_img <- readPNG(png_path)
  img_w <- dim(base_img)[2]
  img_h <- dim(base_img)[1]
  
  barcode_pos <- read.table(gzfile(pos_file), header=FALSE, stringsAsFactors=FALSE)
  colnames(barcode_pos) <- c("Barcode","pos_w","pos_h")
  barcode_pos <- cal_zoom_rate(barcode_pos, img_w, img_h)
  
  df <- data.frame(Barcode=colnames(object),
                   cluster=as.character(object$seurat_clusters)) %>%
    inner_join(barcode_pos, by="Barcode")
  
  cols <- setNames(make_palette(length(unique(df$cluster))), unique(df$cluster))
  
  p_spatial <- ggplot(df, aes(x=pos_w, y=(img_h - pos_h), color=cluster)) +
    background_image(base_img) +
    geom_point(shape=16, size=2.3, alpha=1) +
    coord_fixed(ratio=1, xlim=c(0,img_w), ylim=c(0,img_h), expand=FALSE) +
    scale_color_manual(values=cols) +
    theme_void() +
    labs(color="Cluster", title="Spatial (seurat_clusters)")
  
  ggsave(file.path(out_dir,"spatial_seurat_clusters.png"), p_spatial, width=12, height=10, dpi=300)
  ggsave(file.path(out_dir,"spatial_seurat_clusters.pdf"), p_spatial, width=12, height=10, dpi=600)
} else {
  warning("⚠️ 找不到 HE 图或 barcodes_pos.tsv.gz，跳过空间投影")
}

cat("✅ 绘图完成，结果输出目录：", out_dir, "\n")
