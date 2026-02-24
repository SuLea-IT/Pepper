# ========== 加载库 ==========
library(Seurat)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(png)
library(reshape2)
library(ggpubr)
library(dplyr)
library(tibble)

# 清理环境
rm(list=ls())

# ========== 路径设置 ==========
# 1. Seurat对象RDS文件路径
seurat_rds_path <- "path/to/input_file.rds"

# 2. 空间映射图片路径
png_path = 'path/to/input_file.png'

# 3. barcode位置文件路径
barcode_pos_file = 'path/to/input_file.tsv.gz'

# 4. 输出目录路径
file_dir <- 'path/to/output_dir'

# 5. 基因数据CSV文件路径
gene_data_csv <- 'path/to/input_file.csv'

# 6. 分析标题（用于输出文件夹命名）
analysis_title <- "gene_projection"
# ========================================

# 读取RDS文件
cat("正在读取Seurat对象...\n")
object <- readRDS(seurat_rds_path)
cat("Seurat对象读取成功\n")
cat("- 细胞数:", ncol(object), "\n")
cat("- 基因数:", nrow(object), "\n")

# 确定主要assay
if("SCT" %in% names(object)) {
  primary_assay <- "SCT"
} else if("RNA" %in% names(object)) {
  primary_assay <- "RNA"
} else {
  primary_assay <- names(object)[1]
}
cat("主要使用数据层:", primary_assay, "\n")

# 读取背景图
png <- readPNG(png_path)

# 读取barcode_pos文件
cat("正在读取barcode_pos文件...\n")
barcode_pos_raw <- read.table(gzfile(barcode_pos_file), header = FALSE, stringsAsFactors = FALSE)
barcode_pos <- data.frame(
  Barcode = as.character(barcode_pos_raw[,1]),
  pos_w = as.numeric(barcode_pos_raw[,2]),
  pos_h = as.numeric(barcode_pos_raw[,3]),
  stringsAsFactors = FALSE
)

# ========== 输出目录 ==========
out_path <- paste0(file_dir,'/', analysis_title)
if(!dir.exists(out_path)){dir.create(path=out_path,recursive = TRUE)}
setwd(out_path)

# 读取基因CSV
cat("正在读取基因数据CSV文件...\n")
data_gene <- read.csv(gene_data_csv, stringsAsFactors = FALSE)
cat("基因数据示例:\n")
print(head(data_gene))

# 可用基因列表
if(primary_assay == "SCT"){
  all_available_genes <- rownames(object[["SCT"]])
} else {
  all_available_genes <- rownames(object[["RNA"]])
}

# ========== 基因投影部分 ==========
for(i in 1:nrow(data_gene)){
  gene_id <- data_gene$Gene[i]
  gene_name <- data_gene$GeneName[i]
  cluster <- data_gene$BestCluster[i]
  
  if(!(gene_id %in% all_available_genes)){
    cat("跳过", gene_id, "未找到\n")
    next
  }
  
  if(primary_assay == "SCT"){
    expr_values <- object[["SCT"]]@data[gene_id, ]
  } else {
    expr_values <- object[["RNA"]]@data[gene_id, ]
  }
  
  expr_df <- data.frame(
    Barcode = colnames(object),
    Expr = as.numeric(expr_values),
    stringsAsFactors = FALSE
  )
  
  merged_df <- merge(expr_df, barcode_pos, by = "Barcode", all.x = TRUE)
  
  gene_plot <- ggplot(merged_df, aes(x = pos_w, y = (dim(png)[1] - pos_h))) +
    background_image(png) +
    geom_point(shape = 16, alpha = 1, size = 1.0, aes(color = Expr)) +
    coord_cartesian(xlim = c(0, dim(png)[2]), ylim = c(0, dim(png)[1]), expand = FALSE) +
    scale_color_gradientn(colours = c('#440956', '#299a87', '#fbe725')) +
    theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) +
    ggtitle(paste0(cluster, " - ", gene_name, " (", gene_id, ")"))
  
  # 输出文件名：Cluster_GeneName_Gene.png / .pdf
  out_png <- paste0(out_path, "/", cluster, "_", gene_name, "_", gene_id, ".png")
  out_pdf <- paste0(out_path, "/", cluster, "_", gene_name, "_", gene_id, ".pdf")
  
  ggsave(out_png, plot = gene_plot, width = 10, height = 9, dpi = 300)
  ggsave(out_pdf, plot = gene_plot, width = 10, height = 9, dpi = 300)
  
  cat("完成基因投影:", cluster, gene_name, gene_id, "\n")
}

cat("所有处理完成！\n")
