# ================== 聚类结果空间回投（显示全部 cluster） ==================
library(Seurat)
library(ggplot2)
library(png)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggpubr)   # 用于 background_image

# ========== 需要修改的路径设置 ==========
seurat_rds_path <- "path/to/input_file.rds"
target_sample   <- "25DPA"   # 示例：Fruit_10DPA_ST, Fruit_45DPA_ST, Fruit_25DPA_CA39_ST, ...
png_path        <- "path/to/input_file.png"
barcode_pos_file<- "path/to/input_file.tsv.gz"
file_dir        <- "path/to/output_dir"
analysis_title  <- "cluster_projection"
cluster_column  <- NULL     # 例如: "harmony_clusters", "seurat_clusters", "RNA_snn_res.0.2"

# （可选）cluster 到组织类型的映射（只影响图例标签，不做筛选）
cluster_to_tissue <- list(
)

# ========== 工具函数 ==========
# 去掉样本前缀，例如 "25DPA_L13_0" -> "L13_0"
strip_prefix <- function(x, prefix) {
  gsub(paste0("^", prefix, "_"), "", x)
}

# 计算缩放比例
cal_zoom_rate = function(width, height){
  std_width = 1000
  std_height = std_width / (width * (31)) * (height * (36) * sqrt(3) / 2.0)
  if (std_width / std_height > width / height) {
    scale = width / std_width
  } else {
    scale = height / std_height
  }
  return(scale)
}

# 绘图函数
cluster_spatial_plot <- function(cluster_data, png_img, title = 'Cluster Mapping', point_size = 3.5, cluster_to_tissue = NULL) {
  clv <- sort(unique(as.character(cluster_data$cluster)))
  cluster_data$cluster <- factor(cluster_data$cluster, levels = clv)
  
  custom_colors <- c("0" = "#F56867", "1" = "#FEB915", "2" = "#C798EE")
  colors <- setNames(rep(NA_character_, length(clv)), clv)
  for (nm in intersect(names(custom_colors), clv)) {
    colors[nm] <- custom_colors[nm]
  }
  remaining <- names(colors)[is.na(colors)]
  if (length(remaining) > 0) {
    if (length(remaining) <= 8) {
      default_colors <- RColorBrewer::brewer.pal(max(3, length(remaining)), "Set1")[seq_along(remaining)]
    } else {
      default_colors <- viridis::viridis(length(remaining), option = "D")
    }
    colors[remaining] <- default_colors
  }
  
  # 图例标签
  legend_labels <- setNames(paste0("Cluster ", clv), clv)
  if (!is.null(cluster_to_tissue) && length(cluster_to_tissue) > 0) {
    legend_labels <- setNames(
      ifelse(!is.na(unlist(cluster_to_tissue[clv])),
             paste0("Cluster ", clv, " (", unlist(cluster_to_tissue[clv]), ")"),
             paste0("Cluster ", clv)),
      clv
    )
  }
  # 保证和 colors 对齐
  legend_labels <- legend_labels[names(colors)]
  
  ggplot(cluster_data, aes(x = pos_w, y = (dim(png_img)[1] - pos_h), color = cluster)) +
    background_image(png_img) +
    geom_point(shape = 16, alpha = 0.8, size = point_size) +
    coord_cartesian(xlim = c(0, dim(png_img)[2]), ylim = c(0, dim(png_img)[1]), expand = FALSE) +
    scale_color_manual(values = colors,
                       limits = names(colors),
                       drop = FALSE,
                       labels = legend_labels) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10),
      legend.position = "right"
    ) +
    labs(color = "Tissue Type", title = title)
}

# 简单统计
cluster_statistics <- function(cluster_data) {
  cluster_data %>%
    group_by(cluster) %>%
    summarise(cell_count = n(), .groups = 'drop') %>%
    arrange(cluster)
}

# ========== 主流程 ==========
cat("正在读取Seurat对象...\n")
object <- readRDS(seurat_rds_path)
cat("Seurat对象读取成功\n")

if ("orig.ident" %in% colnames(object@meta.data)) {
  sample_counts <- table(object@meta.data$orig.ident)
  if (!(target_sample %in% names(sample_counts))) {
    stop(paste0("目标样本 ", target_sample, " 不存在。可用样本：", paste(names(sample_counts), collapse = ", ")))
  }
  target_cells <- colnames(object)[object@meta.data$orig.ident == target_sample]
  object_subset <- object[, target_cells]
  cat("目标样本 ", target_sample, " 细胞数：", length(target_cells), "\n")
} else {
  object_subset <- object
  cat("未检测到 orig.ident，使用全部数据\n")
}

# 选择聚类列
meta_columns <- colnames(object@meta.data)
cluster_columns <- grep("cluster", meta_columns, value = TRUE, ignore.case = TRUE)
resolution_columns <- grep("res\\.", meta_columns, value = TRUE, ignore.case = TRUE)
all_cluster_columns <- unique(c(cluster_columns, resolution_columns))

if (!is.null(cluster_column) && cluster_column %in% meta_columns) {
  selected_cluster_col <- cluster_column
} else if ("harmony_clusters" %in% meta_columns) {
  selected_cluster_col <- "harmony_clusters"
} else if ("seurat_clusters" %in% meta_columns) {
  selected_cluster_col <- "seurat_clusters"
} else if (length(all_cluster_columns) > 0) {
  selected_cluster_col <- all_cluster_columns[1]
} else {
  stop("未找到任何聚类信息列！")
}
cat("使用聚类列：", selected_cluster_col, "\n")

# 读 PNG
png_img <- readPNG(png_path)
zoom_scale <- cal_zoom_rate(dim(png_img)[2], dim(png_img)[1])
cat("缩放比例:", zoom_scale, "\n")

# 读 barcode 坐标
barcode_pos_raw <- read.table(gzfile(barcode_pos_file), header = FALSE, stringsAsFactors = FALSE)
barcode_pos <- data.frame(
  Barcode = as.character(barcode_pos_raw[,1]),
  pos_w   = as.numeric(barcode_pos_raw[,2]) * zoom_scale,
  pos_h   = as.numeric(barcode_pos_raw[,3]) * zoom_scale,
  stringsAsFactors = FALSE
)

# 构建聚类数据
cluster_data <- data.frame(
  Barcode = colnames(object_subset),
  cluster = object_subset@meta.data[[selected_cluster_col]],
  stringsAsFactors = FALSE
)

# ========== 条形码对齐（自动去掉前缀） ==========
cluster_data$Barcode_clean <- strip_prefix(cluster_data$Barcode, target_sample)
barcode_pos$Barcode_clean  <- barcode_pos$Barcode

cluster_spatial_data <- merge(
  cluster_data[, c("Barcode", "Barcode_clean", "cluster")],
  barcode_pos[, c("Barcode", "Barcode_clean", "pos_w", "pos_h")],
  by = "Barcode_clean", all.x = TRUE
)

# 保留原始 Seurat 条形码
cluster_spatial_data$Barcode <- cluster_spatial_data$Barcode.x
cluster_spatial_data <- cluster_spatial_data[, c("Barcode", "cluster", "pos_w", "pos_h")]

# 移除无坐标
missing_pos <- sum(is.na(cluster_spatial_data$pos_w) | is.na(cluster_spatial_data$pos_h))
if (missing_pos > 0) {
  cat("警告：有 ", missing_pos, " 个细胞无空间坐标，已移除。\n")
  cluster_spatial_data <- cluster_spatial_data[!is.na(cluster_spatial_data$pos_w) & !is.na(cluster_spatial_data$pos_h), ]
}

cat("最终用于空间映射的细胞数：", nrow(cluster_spatial_data), "\n")

# ========== 输出 ==========
out_path <- file.path(file_dir, analysis_title)
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
setwd(out_path)

main_title <- paste(target_sample, "- Cluster Spatial Mapping (All Clusters)")
cluster_plot <- cluster_spatial_plot(
  cluster_spatial_data, png_img,
  title = main_title, point_size = 6, cluster_to_tissue = cluster_to_tissue
)

output_prefix <- paste0(analysis_title, "_", target_sample)
ggsave(cluster_plot,
       filename = file.path(out_path, paste0(output_prefix, "_cluster_mapping.pdf")),
       width = 14, height = 13)

cat("主要聚类映射图已保存：", file.path(out_path, paste0(output_prefix, "_cluster_mapping.pdf")), "\n")

# 统计表（仅打印）
stats_df <- cluster_statistics(cluster_spatial_data)
stats_df$tissue_type <- sapply(as.character(stats_df$cluster), function(x) {
  if(!is.null(cluster_to_tissue[[x]])) cluster_to_tissue[[x]] else "Unknown"
})
print(stats_df)

cat("\n========== 完成（全量 cluster 显示） ==========\n")
cat("样本: ", target_sample, "\n")
cat("聚类列: ", selected_cluster_col, "\n")
cat("输出目录: ", out_path, "\n")
cat("生成文件: ", paste0(output_prefix, "_cluster_mapping.pdf"), "\n")
