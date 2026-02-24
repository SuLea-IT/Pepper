# ===================== 两个 RDS 合并后重新聚类（自动寻找 parent object） =====================
# v2025-09-01
# - 输入：两个 cluster_*.rds 的路径（或 Seurat 子集 RDS）
# - 自动回溯两级查找各自 object.rds，并按簇号在原对象中截取
# - SCT 集成（SCT Integration）→ PCA/UMAP → 聚类
# - 输出：整合后的 Seurat RDS、UMAP 图（按cluster/样本源）、marker基因等
# ========================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(tibble)
})

## -------------------- 需要你确认的两条路径 --------------------
path_a <- "path/to/input_file.rds"
path_b <- "path/to/input_file.rds"

## -------------------- 整体分析参数 --------------------
res_recluster <- 0.2     # 重新聚类的分辨率
dims_use      <- 1:30    # PCA/UMAP 使用的维度
nfeatures_int <- 3000    # SCT 集成使用的特征数
outdir_base   <- "path/to/output_dir"

## -------------------- 小工具函数 --------------------
get_cluster_id_from_filename <- function(p){
  b <- basename(p)
  if (grepl("cluster_\\d+", b)) {
    sub(".*cluster_(\\d+).*", "\\1", b)
  } else {
    stop("无法从文件名提取 cluster id：", b)
  }
}

get_dpa_from_path <- function(p){
  # 从路径中抓取像 “16DPA”“25DPA” 的片段
  m <- str_extract(p, "(?i)\\b\\d+\\s*DPA\\b")
  if (!is.na(m)) gsub("\\s+", "", toupper(m)) else NA_character_
}

parent_object_path <- function(cluster_rds_path){
  # cluster_*.rds -> .../spatial_projection/ -> 上一级是运行目录 -> 寻找 object.rds
  file.path(dirname(dirname(cluster_rds_path)), "object.rds")
}

# 尝试直接读取 cluster_*.rds；若不是 Seurat 对象，则回溯到 object.rds 再按簇号截取
load_subset_from <- function(cluster_rds_path, add_id = NULL){
  message(">>> 处理: ", cluster_rds_path)
  cl_id <- get_cluster_id_from_filename(cluster_rds_path)
  dpa   <- get_dpa_from_path(cluster_rds_path)
  
  obj_try <- try(readRDS(cluster_rds_path), silent = TRUE)
  if (!inherits(obj_try, "try-error") && inherits(obj_try, "Seurat")){
    message("读取到 Seurat 子集对象：", basename(cluster_rds_path))
    obj_sub <- obj_try
  } else {
    obj_parent_path <- parent_object_path(cluster_rds_path)
    if (!file.exists(obj_parent_path)) {
      stop("未找到父对象 RDS：", obj_parent_path,
           "；请手动把 path_* 指向 Seurat 子集 RDS 或提供父对象路径。")
    }
    message("从父对象截取簇：", obj_parent_path, "  [cluster=", cl_id, "]")
    obj_parent <- readRDS(obj_parent_path)
    
    # 找 cluster 列
    meta_cols <- colnames(obj_parent@meta.data)
    cand <- meta_cols[grepl("cluster", meta_cols, ignore.case = TRUE)]
    cluster_col <- if ("seurat_clusters" %in% meta_cols) "seurat_clusters"
    else if (length(cand) > 0) cand[1]
    else stop("父对象 meta.data 未找到 cluster 列。")
    
    # 截取目标簇
    cl_chr <- as.character(obj_parent@meta.data[[cluster_col]])
    if (!cl_id %in% cl_chr) {
      stop("父对象中未找到目标簇 ", cl_id, "。可用簇：",
           paste(sort(unique(cl_chr)), collapse = ", "))
    }
    obj_sub <- subset(obj_parent, subset = !!as.name(cluster_col) == cl_id)
  }
  
  # 追加源标签
  dpa_lab <- ifelse(is.na(dpa), "Sample", dpa)
  cl_lab  <- paste0("c", get_cluster_id_from_filename(cluster_rds_path))
  src_id  <- if (!is.null(add_id)) add_id else paste0(dpa_lab, "_", cl_lab)
  
  obj_sub$source_batch <- src_id
  # 细胞名加前缀，避免合并后重名
  if (length(grep(paste0("^", src_id, "_"), colnames(obj_sub))) == 0) {
    obj_sub <- RenameCells(obj_sub, add.cell.id = src_id)
  }
  
  list(object = obj_sub, source = src_id)
}

## -------------------- 读取两份子集并构造列表 --------------------
a <- load_subset_from(path_a)  # 返回 list(object, source)
b <- load_subset_from(path_b)

objs <- list(a$object, b$object)
names(objs) <- c(a$source, b$source)
message("子集细胞数：", a$source, "=", ncol(a$object), "；",
        b$source, "=", ncol(b$object))

## -------------------- SCT 集成 --------------------
# 对每个数据集单独做 SCTransform
objs <- lapply(objs, function(x){
  DefaultAssay(x) <- "RNA"
  SCTransform(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = objs, nfeatures = nfeatures_int)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features, verbose = FALSE)

anchors <- FindIntegrationAnchors(object.list = objs,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  verbose = FALSE)

merged <- IntegrateData(anchorset = anchors,
                        normalization.method = "SCT",
                        verbose = FALSE)

DefaultAssay(merged) <- "integrated"

## -------------------- 降维 / 聚类 --------------------
merged <- RunPCA(merged, npcs = max(dims_use), verbose = FALSE)
merged <- RunUMAP(merged, dims = dims_use, verbose = FALSE)
merged <- FindNeighbors(merged, dims = dims_use)
merged <- FindClusters(merged, resolution = res_recluster)

## -------------------- 输出目录 --------------------
tag <- paste0(names(objs), collapse = "_PLUS_")
outdir <- file.path(outdir_base, paste0("merge_", tag, "_SCTint_R", res_recluster))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## -------------------- 可视化 --------------------
# 颜色：cluster 用离散调色板；源样本用另一套
clv <- sort(unique(as.character(merged$seurat_clusters)))
cols_cluster <- setNames(scales::hue_pal()(length(clv)), clv)

p1 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE, repel = TRUE, pt.size = 2.8) + ggtitle("UMAP by cluster")

p2 <- DimPlot(merged, reduction = "umap", group.by = "source_batch",
              pt.size = 2.8) + ggtitle("UMAP by source")

p3 <- DimPlot(merged, reduction = "umap", group.by = "seurat_clusters",
              split.by = "source_batch", ncol = 2, pt.size = 2.8) +
  ggtitle("UMAP by cluster (split by source)")

ggsave(file.path(outdir, "umap_by_cluster.png"), p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "umap_by_source.png"),  p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "umap_split_by_source.png"), p3, width = 12, height = 6, dpi = 300)

## -------------------- 保存对象 & 统计/markers --------------------
saveRDS(merged, file = file.path(outdir, "merged_integrated.rds"))

# 简要统计
stats <- merged@meta.data %>%
  count(source_batch, seurat_clusters, name = "cell_n") %>%
  arrange(source_batch, seurat_clusters)
write.csv(stats, file = file.path(outdir, "cell_counts_by_source_cluster.csv"), row.names = FALSE)

# Marker（整合后 assay 默认 integrated；可切换到 RNA/SCT 再做）
DefaultAssay(merged) <- "integrated"
markers <- FindAllMarkers(merged, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(markers, file = file.path(outdir, "markers_integrated.csv"), row.names = FALSE)

cat("✓ 合并与重新聚类完成。\n输出目录：", outdir, "\n")
