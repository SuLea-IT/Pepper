########################################################################
## Monocle 2 trajectory — 三种 ordering genes 方法对比 (完整 pipeline)
########################################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(ggplot2)
  library(stringr)
  library(Matrix)
})

rm(list = ls())

## 参数
opt <- list(
  outdir = "path/to/output_dir",
  indir  = "path/to/input_file.rds"
)
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE)

## 调色板（针对 cluster_stage3）
base_col <- c("16DPA_all"="#F56867","25_c0"="#FEB915","25_c1"="#C798EE",
              "25_c2"="#59BE86","25_c3"="#7495D3","25_c4"="#D1D1D1")

## 1) 读入 Seurat 对象
obj <- readRDS(opt$indir)

## 如果没有 cluster_stage3，就从 integrated_snn_res.0.3 或 seurat_clusters 生成
if (!"cluster_stage3" %in% colnames(obj@meta.data)) {
  if ("integrated_snn_res.0.3" %in% colnames(obj@meta.data)) {
    obj$cluster_stage3 <- factor(obj$integrated_snn_res.0.3)
  } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    obj$cluster_stage3 <- factor(obj$seurat_clusters)
  } else {
    stop("meta.data 中找不到 cluster_stage3 或常见的分群列，请检查对象")
  }
} else {
  obj$cluster_stage3 <- factor(obj$cluster_stage3)
}

## 2) 提取表达矩阵
assay_use <- if ("SCT" %in% Assays(obj)) "SCT" else "RNA"
expr_matrix <- GetAssayData(obj, assay = assay_use, slot = "counts")

## 3) 构建 phenoData 和 featureData
pd <- new("AnnotatedDataFrame", data = obj@meta.data)
fd <- new("AnnotatedDataFrame",
          data = data.frame(
            gene_short_name = rownames(expr_matrix),
            row.names = rownames(expr_matrix)
          ))

stopifnot(
  nrow(expr_matrix) == nrow(fd),
  all(rownames(expr_matrix) == rownames(fd))
)

## 4) 创建 CellDataSet
mycds0 <- newCellDataSet(
  as(expr_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

## 5) 预处理
mycds0 <- estimateSizeFactors(mycds0)
mycds0 <- estimateDispersions(mycds0)

cat("✅ Seurat 对象已成功转换为 Monocle2 CellDataSet\n")

## ─────────────────── 函数：运行 Monocle 并保存轨迹图 ───────────────── ##
run_monocle <- function(mycds, ordering_genes, method_name, cluster_colors, outdir){
  mycds <- setOrderingFilter(mycds, ordering_genes)
  
  # ordering genes
  p_genes <- plot_ordering_genes(mycds)
  ggsave(file.path(outdir, paste0(method_name,"_genes.png")),
         plot = p_genes, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(method_name,"_genes.pdf")),
         plot = p_genes, width = 6, height = 5, device = cairo_pdf)
  
  mycds <- reduceDimension(mycds, max_components = 2, method = "DDRTree")
  mycds <- orderCells(mycds)
  
  # cluster_stage3
  p1 <- plot_cell_trajectory(mycds, color_by = "cluster_stage3") +
    scale_color_manual(values = cluster_colors)
  ggsave(file.path(outdir, paste0(method_name,"_cluster.png")),
         plot = p1, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(method_name,"_cluster.pdf")),
         plot = p1, width = 6, height = 5, device = cairo_pdf)
  
  # Pseudotime
  p2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
  ggsave(file.path(outdir, paste0(method_name,"_pseudotime.png")),
         plot = p2, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(method_name,"_pseudotime.pdf")),
         plot = p2, width = 6, height = 5, device = cairo_pdf)
  
  # State
  p3 <- plot_cell_trajectory(mycds, color_by = "State")
  ggsave(file.path(outdir, paste0(method_name,"_state.png")),
         plot = p3, width = 6, height = 5, dpi = 300)
  ggsave(file.path(outdir, paste0(method_name,"_state.pdf")),
         plot = p3, width = 6, height = 5, device = cairo_pdf)
  
  return(mycds)
}

## ─────────────────── 方法 1: Monocle 高离散度基因 ───────────────── ##
disp_table <- dispersionTable(mycds0)
disp_genes <- subset(disp_table,
                     mean_expression >= 0.1 &
                       dispersion_empirical >= 1 * dispersion_fit)$gene_id

mycds_disp <- run_monocle(mycds0, disp_genes, "monocle_dispersion",
                          cluster_colors = base_col, outdir = opt$outdir)

## ─────────────────── 方法 2: Seurat 高变基因 ───────────────────── ##
var_genes <- VariableFeatures(obj)
mycds_var <- run_monocle(mycds0, var_genes, "seurat_variable",
                         cluster_colors = base_col, outdir = opt$outdir)

cat("🎉 三种 ordering genes 方法完成；结果保存在：", opt$outdir, "\n")

## ─────────────────── 8. 导出指定基因的拟时序曲线 ───────────────── ##
genes_to_plot <- c("LOC107859694" = "CS",
                   "LOC107878306" = "MYB31",
                   "LOC107879881" = "WRKY9",
                   "LOC107864105" = "Kasl",
                   "LOC107875336" = "COMT36")

min_expr_tbl <- c(CS = 0, MYB31 = 0, WRKY9 = 0, Kasl = 0, COMT36 = 0)

mycds_final <- mycds_var  # 选择一个结果用于绘图

# ========== A. 普通拟时序（所有分支合并） ==========
for (gene_id in names(genes_to_plot)) {
  gene_sym <- genes_to_plot[gene_id]
  min_e    <- min_expr_tbl[gene_sym]
  
  p_gene <- plot_genes_in_pseudotime(mycds_final[gene_id, ],
                                     color_by  = "cluster_stage3",
                                     min_expr  = min_e,
                                     cell_size = 2) +
    scale_color_manual(values = base_col) +
    ggtitle(gene_sym)
  
  ggsave(file.path(opt$outdir, paste0(gene_sym, "_pseudotime.png")),
         plot = p_gene, width = 6, height = 4.5, dpi = 300)
  ggsave(file.path(opt$outdir, paste0(gene_sym, "_pseudotime.pdf")),
         plot = p_gene, width = 6, height = 4.5, device = cairo_pdf)
}

# ========== B. 两个分支对比（实线/虚线） ==========
plot_gene_two_directions <- function(cds, gene_id, gene_sym,
                                     states = c(2, 3),
                                     min_expr = 0.0,
                                     cluster_cols = NULL,
                                     point_size = 2,
                                     span = 0.6) {
  pd <- pData(cds)
  stopifnot(all(c("Pseudotime","State","cluster_stage3") %in% colnames(pd)))
  
  expr_v <- log2(as.numeric(exprs(cds[gene_id, ])) + 1)
  
  df <- data.frame(
    Pseudotime = pd$Pseudotime,
    Expr       = expr_v,
    State      = pd$State,
    Cluster    = pd$cluster_stage3,
    stringsAsFactors = FALSE
  )
  
  df <- df[df$State %in% states & !is.na(df$Pseudotime) & df$Expr >= min_expr, , drop = FALSE]
  df$Direction <- factor(df$State, levels = states,
                         labels = paste0("State ", states))
  
  p <- ggplot(df, aes(Pseudotime, Expr)) +
    geom_point(aes(color = Cluster), size = point_size, alpha = 0.65) +
    geom_smooth(aes(linetype = Direction),
                method = "loess", se = FALSE, span = span,
                color = "black", size = 1.1) +
    scale_color_manual(values = cluster_cols) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(title = paste0(gene_sym, " (branched pseudotime)"),
         x = "Pseudotime", y = "Expression",
         linetype = "Direction", color = "Cluster") +
    theme_classic(base_size = 12)
  
  return(p)
}

for (gene_id in names(genes_to_plot)) {
  gene_sym <- genes_to_plot[gene_id]
  min_e    <- min_expr_tbl[gene_sym]
  
  p_gene <- plot_gene_two_directions(
    cds          = mycds_final,
    gene_id      = gene_id,
    gene_sym     = gene_sym,
    states       = c(2, 3),    # ← 修改为你实际想比较的 State 编号
    min_expr     = min_e,
    cluster_cols = base_col,
    point_size   = 2,
    span         = 0.6
  )
  
  ggsave(file.path(opt$outdir, paste0(gene_sym, "_branched.png")),
         plot = p_gene, width = 6, height = 4.5, dpi = 300)
  ggsave(file.path(opt$outdir, paste0(gene_sym, "_branched.pdf")),
         plot = p_gene, width = 6, height = 4.5, device = cairo_pdf)
}

cat("✅ 基因拟时序曲线 (普通版 + Branched版) 导出完成；结果保存在：", opt$outdir, "\n")
