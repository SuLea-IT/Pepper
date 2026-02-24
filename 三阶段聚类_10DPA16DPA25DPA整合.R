# ========== 三阶段聚类分析: 10DPA整体 vs 16DPA整体 vs 25DPA子簇 ==========
# v2025-09-02.stage123_SCT_selectfeatures
suppressPackageStartupMessages({
  library(Seurat); library(sctransform)
  library(dplyr); library(ggplot2); library(tibble); library(scales)
})

set.seed(1234)

# ---------------- 参数 ----------------
rds_10DPA <- "path/to/input_file.rds"
rds_16DPA <- "path/to/input_file.rds"
rds_25DPA <- "path/to/input_file.rds"
res2      <- 0.2     # Stage2 resolution
npca      <- 30
outdir <- "path/to/output_dir"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---------------- 1) 读取三个 RDS ----------------
obj10 <- readRDS(rds_10DPA); obj10$sample <- "10DPA"
obj16 <- readRDS(rds_16DPA); obj16$sample <- "16DPA"
obj25 <- readRDS(rds_25DPA); obj25$sample <- "25DPA"

# ---------------- 确保三个对象都有 SCT assay ----------------
if (!"SCT" %in% Assays(obj10)) obj10 <- SCTransform(obj10, verbose = FALSE)
if (!"SCT" %in% Assays(obj16)) obj16 <- SCTransform(obj16, verbose = FALSE)
if (!"SCT" %in% Assays(obj25)) obj25 <- SCTransform(obj25, verbose = FALSE)

# ---------------- Stage1: 样本分簇 ----------------
obj.all <- merge(obj10, y = list(obj16, obj25))
DefaultAssay(obj.all) <- "SCT"

# 使用 SelectIntegrationFeatures 选基因
features <- SelectIntegrationFeatures(object.list = list(obj10, obj16, obj25), nfeatures = 3000)
VariableFeatures(obj.all) <- features

# PCA + UMAP
obj.all <- RunPCA(obj.all, npcs = npca, verbose = FALSE)
obj.all <- RunUMAP(obj.all, dims = 1:npca)

obj.all$cluster_stage1 <- obj.all$sample

p1 <- DimPlot(obj.all, group.by = "cluster_stage1", label = TRUE, repel = TRUE) +
  ggtitle("Stage1: 10DPA vs 16DPA vs 25DPA")
ggsave(file.path(outdir, "stage1_umap.png"), p1, width = 6, height = 5, dpi = 300)

# ---------------- Stage2: 25DPA 内部分簇 ----------------
obj25.sub <- subset(obj.all, sample == "25DPA")
DefaultAssay(obj25.sub) <- "SCT"

# 这里重新挑 variable features
features.sub <- SelectIntegrationFeatures(object.list = list(obj25.sub), nfeatures = 3000)
VariableFeatures(obj25.sub) <- features.sub

obj25.sub <- RunPCA(obj25.sub, npcs = npca, verbose = FALSE)
obj25.sub <- FindNeighbors(obj25.sub, dims = 1:npca)
obj25.sub <- FindClusters(obj25.sub, resolution = res2)
obj25.sub <- RunUMAP(obj25.sub, dims = 1:npca)

obj25.sub$cluster_stage2 <- paste0("25_c", obj25.sub$seurat_clusters)

p2 <- DimPlot(obj25.sub, group.by = "cluster_stage2", label = TRUE, repel = TRUE) +
  ggtitle("Stage2: 25DPA subclusters")
ggsave(file.path(outdir, "stage2_umap.png"), p2, width = 6, height = 5, dpi = 300)

# ---------------- Stage3: 合并 10DPA整体 + 16DPA整体 + 25DPA子簇 ----------------
obj10$cluster_stage3 <- "10DPA_all"
obj16$cluster_stage3 <- "16DPA_all"
obj25.sub$cluster_stage3 <- obj25.sub$cluster_stage2

obj.final <- merge(obj10, y = list(obj16, obj25.sub))
DefaultAssay(obj.final) <- "SCT"

# 重新选 features
features.final <- SelectIntegrationFeatures(object.list = list(obj10, obj16, obj25.sub), nfeatures = 3000)
VariableFeatures(obj.final) <- features.final

obj.final <- RunPCA(obj.final, npcs = npca, verbose = FALSE)
obj.final <- FindNeighbors(obj.final, dims = 1:npca)
obj.final <- FindClusters(obj.final, resolution = 0)  # res=0, 不再分裂
obj.final <- RunUMAP(obj.final, dims = 1:npca)

p3 <- DimPlot(obj.final, group.by = "cluster_stage3", label = TRUE, repel = TRUE) +
  ggtitle("Stage3: 10DPA_all vs 16DPA_all vs 25DPA subclusters")
ggsave(file.path(outdir, "stage3_umap.png"), p3, width = 6, height = 5, dpi = 300)

# ---------------- 保存对象 ----------------
saveRDS(obj.all,   file.path(outdir, "stage1_obj_all.rds"))
saveRDS(obj25.sub, file.path(outdir, "stage2_obj25_sub.rds"))
saveRDS(obj.final, file.path(outdir, "stage3_obj_final.rds"))

cat("✅ 完成三阶段聚类分析 (SCT + SelectIntegrationFeatures)；结果目录：", outdir, "\n")
