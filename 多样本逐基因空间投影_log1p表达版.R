# ===================== 逐基因投影（log1p(counts)，0=固定色，正值渐变；对齐 cal_zoom_rate + 指定点大小） =====================

library(Seurat)
library(tidyverse)
library(png)
library(ggplot2)
library(ggpubr)   # background_image
library(Matrix)

## ---------- 配置 ----------
gene_set_csv_path   <- "path/to/input_file.csv"  # 需包含 gene_id,gene_name 两列
base_dir            <- "path/to/input_dir"
output_base_dir     <- "path/to/output_dir"
folders_to_process  <- c("Fruit_10DPA","Fruit_16DPA","Fruit_25DPA","Fruit_35DPA","Fruit_45DPA")
dpa_level_mapping   <- list("Fruit_10DPA"="L7","Fruit_16DPA"="L13","Fruit_25DPA"="L13","Fruit_35DPA"="L9","Fruit_45DPA"="L13")
data_subdir_template<- "Data/subdata/{LEVEL}_heAuto/"
png_subpath         <- "Data/he_roi_small.png"
barcode_pos_filename<- "barcodes_pos.tsv.gz"

# 若没有“直接数据”文件夹，保持空向量
direct_data_folders <- character(0)

# 默认绘图参数（未知 Level 的兜底值）
point_size_default <- 1.0
point_alpha        <- 1
point_shape        <- 16
image_dpi          <- 300

# 三色设置
col_zero <- "#99c7de"  # 0 表达固定色
col_low  <- "#f2f222"  # 低表达
col_high <- "#ca4f80"  # 高表达

# 可选微调：若仍有轻微偏移，可设置像素偏移（正值向右/向上）
x_offset <- 0
y_offset <- 0

dir.create(output_base_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- 读取基因表（两列：gene_id, gene_name） ----------
gene_df <- read.csv(gene_set_csv_path, stringsAsFactors = FALSE)
if (!all(c("gene_id","gene_name") %in% colnames(gene_df))) {
  stop("CSV 必须包含列：gene_id, gene_name")
}
cat("输入基因条目数：", nrow(gene_df), "\n")

## ---------- 缩放函数（与参考流程一致） ----------
cal_zoom_rate <- function(width, height){
  std_width  <- 1000
  std_height <- std_width / (46 * 31) * (46 * 36 * sqrt(3) / 2.0)
  if (std_width / std_height > width / height) width / std_width else height / std_height
}

## ---------- 工具函数 ----------
get_level_for_dpa <- function(dpa_name) {
  if (dpa_name %in% names(dpa_level_mapping)) dpa_level_mapping[[dpa_name]] else "L13"
}
get_data_subdir <- function(folder_name) {
  if (folder_name %in% direct_data_folders) "" else gsub("\\{LEVEL\\}", get_level_for_dpa(folder_name), data_subdir_template)
}
# 条码自适配：第一轮原样匹配；若全 NA，则尝试去掉或补上 `-1` 尾缀
smart_join_expr <- function(spatial_coords, expr_df) {
  joined <- left_join(spatial_coords, expr_df, by = "Barcode")
  if (all(is.na(joined$expr))) {
    if (any(endsWith(expr_df$Barcode, "-1"))) {
      expr_df2 <- mutate(expr_df, Barcode = sub("-1$", "", Barcode))
    } else {
      expr_df2 <- mutate(expr_df, Barcode = paste0(Barcode, "-1"))
    }
    joined <- left_join(spatial_coords, expr_df2, by = "Barcode")
  }
  joined
}
# 安全文件名
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)

# —— 固定点大小映射：不再计算 ——
level_pointsize <- c("L7" = 0.8, "L9" = 0.7, "L13" = 1.2)

## ---------- 主循环（逐样本 × 逐基因） ----------
for (folder in folders_to_process) {
  cat(">>> 样本：", folder, "\n")
  
  # 当前 Level 与点大小（直接查表，未知则用兜底值）
  level_str  <- get_level_for_dpa(folder)
  pt_size    <- if (level_str %in% names(level_pointsize)) level_pointsize[[level_str]] else point_size_default
  cat("   - Level:", level_str, "-> point_size:", pt_size, "\n")
  
  data_dir <- file.path(base_dir, folder, get_data_subdir(folder))
  png_path <- if (folder %in% direct_data_folders) file.path(base_dir, folder, "he_roi_small.png")
  else file.path(base_dir, folder, png_subpath)
  barcode_pos_path <- file.path(data_dir, barcode_pos_filename)
  
  # 输出子目录（每个样本一个子目录）
  out_dir <- file.path(output_base_dir, folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 表达矩阵（兼容 10X 返回 list）
  data10x <- try(Read10X(data.dir = data_dir), silent = TRUE)
  if (inherits(data10x, "try-error")) { cat("   × 跳过，Read10X 失败：", data_dir, "\n"); next }
  if (is.list(data10x)) {
    data10x <- if ("Gene Expression" %in% names(data10x)) data10x[["Gene Expression"]] else data10x[[1]]
  }
  obj <- CreateSeuratObject(counts = data10x)
  
  # HE 图像 & 条码坐标
  png_img <- try(readPNG(png_path), silent = TRUE)
  if (inherits(png_img, "try-error")) { cat("   × 跳过，PNG 不存在：", png_path, "\n"); next }
  img_w <- dim(png_img)[2]; img_h <- dim(png_img)[1]
  
  barcode_pos <- try(
    read.table(gzfile(barcode_pos_path), header = FALSE) %>%
      dplyr::rename(Barcode = V1, pos_w = V2, pos_h = V3),
    silent = TRUE
  )
  if (inherits(barcode_pos, "try-error")) { cat("   × 跳过，坐标不存在：", barcode_pos_path, "\n"); next }
  
  # 与参考流程一致的缩放
  zoom_scale <- cal_zoom_rate(img_w, img_h)
  spatial_coords <- barcode_pos %>%
    mutate(
      pos_w = pos_w * zoom_scale + x_offset,
      pos_h = pos_h * zoom_scale + y_offset
    )
  
  # 预取 counts 矩阵与行名
  mat_counts <- GetAssayData(obj, slot = "counts")
  rn <- rownames(mat_counts)
  
  # 逐基因绘图
  for (i in seq_len(nrow(gene_df))) {
    gid <- gene_df$gene_id[i]
    gnm <- gene_df$gene_name[i]
    
    # 优先使用 gene_id，其次 gene_name
    row_use <- NA_character_
    if (!is.na(gid) && gid %in% rn) row_use <- gid else if (!is.na(gnm) && gnm %in% rn) row_use <- gnm
    if (is.na(row_use)) {
      cat("   - 跳过：", gid, "/", gnm, "（未在表达矩阵中找到）\n")
      next
    }
    
    # 表达向量（>0 渐变强度用 log1p）
    expr <- as.numeric(mat_counts[row_use, ])
    names(expr) <- colnames(mat_counts)
    expr_df <- data.frame(Barcode = names(expr), expr = expr, stringsAsFactors = FALSE)
    
    # 合并（含 -1 自适配）
    plot_df <- smart_join_expr(spatial_coords, expr_df)
    if (!nrow(plot_df)) { cat("   - 跳过：条码不匹配\n"); next }
    
    # 分层（0 固定色；>0 渐变）
    plot_df0 <- dplyr::filter(plot_df, expr == 0)
    plot_df1 <- dplyr::filter(plot_df, expr  > 0)
    if (nrow(plot_df0) == 0 && nrow(plot_df1) == 0) { cat("   - 跳过：无数据\n"); next }
    if (nrow(plot_df1) > 0) plot_df1$expr_pos <- log1p(plot_df1$expr)
    
    gene_label <- ifelse(!is.na(gnm) && nzchar(gnm), paste0(gnm, " (", gid, ")"), gid)
    file_stub  <- safe_name(paste0(gid, "_", ifelse(is.na(gnm) || !nzchar(gnm), "NA", gnm)))
    out_png    <- file.path(out_dir, paste0(folder, "_", file_stub, "_基因投影.png"))
    
    p <- ggplot() +
      background_image(png_img) +
      { if (nrow(plot_df0) > 0)
        geom_point(
          data = plot_df0,
          aes(x = pos_w, y = (img_h - pos_h)),
          color = col_zero, shape = point_shape, size = pt_size, alpha = point_alpha
        )
        else NULL } +
      { if (nrow(plot_df1) > 0)
        geom_point(
          data = plot_df1,
          aes(x = pos_w, y = (img_h - pos_h), color = expr_pos),
          shape = point_shape, size = pt_size, alpha = point_alpha
        )
        else NULL } +
      scale_color_gradient(low = col_low, high = col_high, name = "表达\nlog1p") +
      coord_cartesian(xlim = c(0, img_w), ylim = c(0, img_h), expand = FALSE) +
      theme_void() +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(folder, " - ", gene_label))
    
    ggsave(out_png, plot = p, width = img_w + 200, height = img_h, units = "px", dpi = image_dpi)
    cat("   ✓ 已保存：", out_png, "\n")
  }
}

cat("全部完成。\n")
