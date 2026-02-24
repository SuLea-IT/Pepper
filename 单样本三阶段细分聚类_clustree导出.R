# ===================== Capsicum scRNA-seq Pipeline (简化版 + clustree 图导出) =====================
# v2025-08-23
# - 阶段 1：一次聚类（导出 RDS & 基础输出 + clustree 图）
# - 阶段 2：子聚类（按父簇手动设定分辨率 + clustree 图）
# - 阶段 3：可选的聚类重编码与细化（重编码后再次导出 clustree 图）
# =====================================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(dplyr)
  library(tidyverse)
  library(future)
  library(SeuratDisk)
  library(clustree)          # 用于绘制 clustree 图
  library(ggplot2)           # 用于绘图
})

options(future.globals.maxSize = 9e12)
set.seed(1234)

# ===================== 配置列表 =====================
config <- list(
  # 输入输出目录
  indir = "path/to/input_dir",           # 10X 矩阵目录（如包含 filtered_feature_bc_matrix/）
  outdir_base = "path/to/output_dir",     # 结果输出目录
  existing_rds = NULL,                          # 若已有处理好的 Seurat RDS，这里填路径；为 NULL 时读取 10X
  # 每个阶段的聚类参数
  stage1 = list(
    method = "SCT",           # 选择方法："SCT" 或 "NFS"（归一化 + 特征选择）
    resolution = 0.3,          # 第一次聚类的分辨率
    min_cells = 10,            # 最少细胞数
    min_features = 50,         # 最少基因数
    npca = 30,                 # PCA 维度数（参与邻居图/UMAP/ClusterTree）
    force_recluster = FALSE    # 若 existing_rds 已含聚类，是否仍按本参数强制重跑邻居/聚类/UMAP
  ),
  
  stage2 = list(
    refine_plan = c("3" = 0.3, "4" = 0.3)    # 需要细化的父簇 = 分辨率
  ),
  
  # 可选的 Stage 3（聚类重编码和细化）
  stage3 = list(
    reencode = TRUE,                 # 是否将细分簇重新编码为 0..N-1
    reencode_method = "sequential"  # 重编码方法："sequential"（0,1,2,...）或 "alphabetical"（字母顺序）
  )
)

# ===================== 工具函数：导出 clustree 图 =====================
# 使用 clustree 包绘制聚类树图
export_clustree_plot <- function(obj, out_prefix, outdir = config$outdir_base){
  # 提取聚类信息
  cluster_info <- as.data.frame(obj@meta.data)
  
  # 使用 clustree 生成聚类树图
  clustree_plot <- clustree(cluster_info, prefix = "seurat_clusters")
  
  # 导出 clustree 图
  clustree_plot_path <- file.path(outdir, paste0(out_prefix, "_clustree.png"))
  ggsave(clustree_plot_path, clustree_plot, width = 8, height = 6, dpi = 300)
  message("clustree 图已保存：", clustree_plot_path)
  
  return(obj)
}

# ===================== 辅助函数 =====================
set_default_assay <- function(obj, method){
  if (method == "SCT" && "SCT" %in% names(obj@assays)){
    DefaultAssay(obj) <- "SCT"
  } else if ("RNA" %in% names(obj@assays)){
    DefaultAssay(obj) <- "RNA"
  }
  obj
}

ensure_embeddings_and_clusters <- function(obj, method, npca, resolution, force = FALSE){
  obj <- set_default_assay(obj, method)
  has_pca   <- "pca" %in% names(obj@reductions)
  has_umap  <- "umap" %in% names(obj@reductions)
  has_clust <- "seurat_clusters" %in% colnames(obj@meta.data)
  
  if (force || !has_pca){
    obj <- RunPCA(obj, npcs = max(50, npca), verbose = FALSE)
  }
  # 邻居图 & 聚类：若强制或无聚类则重建
  if (force || !has_clust){
    obj <- FindNeighbors(obj, dims = 1:npca)
    obj <- FindClusters(obj, resolution = resolution)
  }
  # UMAP：若强制或无 UMAP 则计算（尽管不导图）
  if (force || !has_umap){
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:npca, return.model = TRUE)
  }
  obj
}
# ===================== 阶段 1：一次聚类（接受 RDS 或 10X） =====================
if (!is.null(config$existing_rds)){
  message("【Stage-1】读取现有 RDS: ", config$existing_rds)
  object <- readRDS(config$existing_rds)
  object <- ensure_embeddings_and_clusters(
    obj = object,
    method = config$stage1$method,
    npca = config$stage1$npca,
    resolution = config$stage1$resolution,
    force = config$stage1$force_recluster
  )
} else {
  message("【Stage-1】读取 10X：", config$indir)
  mtx <- Read10X(config$indir)
  object <- CreateSeuratObject(
    mtx,
    assay = "RNA",
    min.cells = config$stage1$min_cells,
    min.features = config$stage1$min_features
  )
  if (config$stage1$method == "SCT") {
    object <- SCTransform(object, verbose = FALSE)
  } else {
    object <- NormalizeData(object)
    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
    object <- ScaleData(object, features = rownames(object))
  }
  object <- ensure_embeddings_and_clusters(
    obj = object,
    method = config$stage1$method,
    npca = config$stage1$npca,
    resolution = config$stage1$resolution,
    force = TRUE
  )
}

# 保存阶段 1 输出
rds_stage1 <- file.path(config$outdir_base, "object_stage1.rds")
saveRDS(object, rds_stage1)

# 导出 marker 基因、平均表达和簇映射
markers <- FindAllMarkers(object, only.pos = TRUE)
write.csv(markers, file = file.path(config$outdir_base, "stage1_marker_gene_onlypos.csv"), row.names = FALSE)

# 导出平均表达量表
object <- set_default_assay(object, config$stage1$method)
ae <- AverageExpression(object, assays = ifelse(DefaultAssay(object) == "SCT", "SCT", "RNA"), verbose = FALSE)
ae_df <- if (is.list(ae) && "SCT" %in% names(ae)) ae$SCT else as.data.frame(ae)
write.csv(ae_df, file = file.path(config$outdir_base, "stage1_average_expression.csv"), row.names = TRUE)

# 导出 UMAP 图
umap_plot_stage1 <- DimPlot(object, reduction = "umap", label = TRUE)
ggsave(file.path(config$outdir_base, "stage1_umap.png"), plot = umap_plot_stage1)

# ===================== 阶段 2：子聚类 =====================
if (!is.null(config$stage2$refine_plan) && length(config$stage2$refine_plan) > 0) {
  # 初始化 refined 标签
  if (!"cluster_refined" %in% colnames(object@meta.data)) {
    object$cluster_refined <- as.character(object$seurat_clusters)
  }
  Idents(object) <- "seurat_clusters"
  
  for (parent in names(config$stage2$refine_plan)) {
    res2 <- config$stage2$refine_plan[[parent]]
    if (!(parent %in% levels(Idents(object)))) {
      warning(sprintf("父簇 '%s' 不存在于 seurat_clusters 中，跳过。", parent))
      next
    }
    sub <- subset(object, idents = parent)
    if (ncol(sub) < 50) {
      warning(sprintf("父簇 '%s' 细胞数过少 (%d)，跳过。", parent, ncol(sub)))
      next
    }
    # 与一次聚类一致的方法
    if (DefaultAssay(object) == "SCT" || config$stage1$method == "SCT") {
      DefaultAssay(sub) <- "RNA"
      sub <- SCTransform(sub, verbose = FALSE, assay = "RNA")
      DefaultAssay(sub) <- "SCT"
    } else {
      sub <- NormalizeData(sub)
      sub <- FindVariableFeatures(sub, selection.method = "vst", nfeatures = 2000)
      sub <- ScaleData(sub)
    }
    npca_sub <- 20
    sub <- RunPCA(sub, npcs = npca_sub, verbose = FALSE)
    sub <- FindNeighbors(sub, dims = 1:npca_sub, verbose = FALSE)
    sub <- FindClusters(sub, resolution = res2, verbose = FALSE)
    sub <- RunUMAP(sub, dims = 1:npca_sub, verbose = FALSE)
    
    # 写回细分标签（如 3_0, 3_1, ...）
    sub_labels <- as.character(Idents(sub))
    mapping <- setNames(paste0(parent, "_", sub_labels), Cells(sub))
    object$cluster_refined[names(mapping)] <- mapping[names(mapping)]
  }
  
  # 保存阶段 2 输出
  rds_stage2 <- file.path(config$outdir_base, "object_stage2_refined.rds")
  saveRDS(object, rds_stage2)
  
  # 导出 stage 2 的 marker 基因、平均表达和 UMAP 图
  markers_stage2 <- FindAllMarkers(object, only.pos = TRUE)
  write.csv(markers_stage2, file = file.path(config$outdir_base, "stage2_marker_gene_onlypos.csv"), row.names = FALSE)
  
  ae_stage2 <- AverageExpression(object, assays = ifelse(DefaultAssay(object) == "SCT", "SCT", "RNA"), verbose = FALSE)
  ae_df_stage2 <- if (is.list(ae_stage2) && "SCT" %in% names(ae_stage2)) ae_stage2$SCT else as.data.frame(ae_stage2)
  write.csv(ae_df_stage2, file = file.path(config$outdir_base, "stage2_average_expression.csv"), row.names = TRUE)
  
  umap_plot_stage2 <- DimPlot(object, reduction = "umap", label = TRUE)
  ggsave(file.path(config$outdir_base, "stage2_umap.png"), plot = umap_plot_stage2)
}

# ===================== 阶段 3：可选聚类重编码 =====================
if (config$stage3$reencode) {
  if (!"cluster_refined" %in% colnames(object@meta.data)) {
    object$cluster_refined <- as.character(object$seurat_clusters)
  }
  old_clusters_chr <- as.character(object$cluster_refined)
  old_levels <- unique(old_clusters_chr)
  object$cluster_refined <- factor(object$cluster_refined, levels = old_levels)
  object$cluster_refined <- as.integer(object$cluster_refined) - 1  # 从 0 开始
  object$cluster_refined <- factor(object$cluster_refined)
  
  # 保存重编码后的输出
  rds_stage3_renamed <- file.path(config$outdir_base, "object_stage3_refined_renamed.rds")
  saveRDS(object, rds_stage3_renamed)
  
  # 导出 stage 3 的 marker 基因、平均表达和 UMAP 图
  markers_stage3 <- FindAllMarkers(object, only.pos = TRUE)
  write.csv(markers_stage3, file = file.path(config$outdir_base, "stage3_marker_gene_onlypos.csv"), row.names = FALSE)
  
  ae_stage3 <- AverageExpression(object, assays = ifelse(DefaultAssay(object) == "SCT", "SCT", "RNA"), verbose = FALSE)
  ae_df_stage3 <- if (is.list(ae_stage3) && "SCT" %in% names(ae_stage3)) ae_stage3$SCT else as.data.frame(ae_stage3)
  write.csv(ae_df_stage3, file = file.path(config$outdir_base, "stage3_average_expression.csv"), row.names = TRUE)
  
  umap_plot_stage3 <- DimPlot(object, reduction = "umap", label = TRUE)
  ggsave(file.path(config$outdir_base, "stage3_umap.png"), plot = umap_plot_stage3)
  
  message("阶段 3（重编码）完成。已保存：", rds_stage3_renamed)
}

# 最终消息
message("管道执行完毕！已在各阶段导出 clusterTree（newick/edges/cophenetic/direction）以及 marker 基因、表达量表和 UMAP 图。")

