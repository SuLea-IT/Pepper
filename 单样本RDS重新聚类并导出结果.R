# 20240429更新，添加了单样品的loom文件生成
# s1.2_20240513更新，添加不同聚类cluster导出
# 20250901更新，改为从 RDS 读取 Seurat 对象

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(SeuratDisk)
  library(dplyr)
  library(tibble)
})

set.seed(1234)

###### 参数设置
opt <- list(
  # res为resolution,值越大聚类个数越多
  res = 0.3,
  # 以下两个参数在读取RDS时不用于过滤，仅保留记录
  MinCell = 10,
  MinFeatures = 50,
  # PCA维度数量
  npca = 30,
  # 数据标准化的方法 "SCT" 或 "NFS"
  method = "SCT",
  # 读取已有的 Seurat RDS 对象（请修改为你的RDS路径）
  in_rds = "path/to/input_file.rds",
  # 结果输出路径（会自动在内部拼接方法与参数）
  outdir = "path/to/output_dir"
)

# 自定义每个cluster的颜色
col <- c(
  "0"="#F56867","1"="#FEB915","2"="#C798EE","3"="#59BE86","4"="#7495D3","5"="#D1D1D1",
  "6"="#6D1A9C","7"="#15821E","8"="#3A84E6","9"="#70e014","10"="#787878","11"="#DB4C6C",
  "12"="#0430e0","13"="#554236","14"="#AF5F3C","15"="#ff7700","16"="#e00417","17"="#DAB370",
  "18"="#fcfc05","19"="#268785","20"="#09f9f5","21"="#246b93","22"="#cc8e12","23"="#d561dd",
  "24"="#c93f00","25"="#ddd53e","26"="#4aef7b","27"="#e86502","28"="#9ed84e","29"="#39ba30",
  "30"="#6ad157","31"="#8249aa","32"="#99db27","33"="#e07233","34"="#ff523f"
)

#### 参数传递
res         <- opt$res
MinCell     <- opt$MinCell
MinFeatures <- opt$MinFeatures
npca        <- opt$npca
method      <- opt$method
in_rds      <- opt$in_rds
outdir      <- paste0(opt$outdir, "/", method, "_C", MinCell, "_F", MinFeatures, "_P", npca, "_R", res)
if (!dir.exists(outdir)) dir.create(path = outdir, recursive = TRUE)

#### 从RDS读取 Seurat 对象
cat("Reading Seurat object from RDS...\n")
object <- readRDS(in_rds)
stopifnot(inherits(object, "Seurat"))

#### 数据标准化 / 变量基因 / 缩放
if (method == "SCT") {
  message("Method: SCT")
  if (!"SCT" %in% Assays(object)) {
    # 若尚无SCT assay，则从RNA或默认assay进行SCTransform
    base_assay <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
    DefaultAssay(object) <- base_assay
    object <- SCTransform(object, verbose = FALSE, assay = base_assay)
  }
  DefaultAssay(object) <- "SCT"
} else {
  message("Method: NFS (Normalize-Feature-Scale)")
  # 以 RNA 为默认assay（若存在）
  if ("RNA" %in% Assays(object)) DefaultAssay(object) <- "RNA"
  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  object <- ScaleData(object, features = rownames(object), verbose = FALSE)
}

#### PCA降维
object <- RunPCA(object, npcs = 50, verbose = FALSE, reduction.name = "pca")
plot_elb <- ElbowPlot(object, ndims = 50)
ggsave(filename = file.path(outdir, "elbowplot.png"), plot = plot_elb, width = 10, height = 8)

#### 细胞聚类与UMAP
object <- FindNeighbors(object, dims = 1:npca, verbose = FALSE)
object <- FindClusters(object, resolution = res, verbose = FALSE)
object <- RunUMAP(object, reduction = "pca", dims = 1:npca, verbose = FALSE)

umapplot <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters", pt.size = 2, cols = col)
ggsave(filename = file.path(outdir, "umap_cluster.png"), plot = umapplot, width = 10, height = 8)

#### 发现特征基因（only.pos=TRUE）
object.markers <- FindAllMarkers(object, only.pos = TRUE, verbose = FALSE)
write.csv(object.markers, file = file.path(outdir, "marker_gene_onlypos.csv"), row.names = FALSE)

#### 拟bulk基因平均表达
if (method == "SCT" && "SCT" %in% Assays(object)) {
  avg_mat <- AverageExpression(object, assays = "SCT", slot = "data", verbose = FALSE)$SCT
} else {
  # 回退到 RNA
  assay_for_avg <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  avg_mat <- AverageExpression(object, assays = assay_for_avg, slot = "data", verbose = FALSE)[[assay_for_avg]]
}
write.csv(avg_mat, file = file.path(outdir, "average_expression.csv"))

#### 保存RDS（处理后的对象）
saveRDS(object, file = file.path(outdir, "object.rds"))


#### 导出聚类信息
Cluster <- object@meta.data %>%
  dplyr::select(seurat_clusters) %>%
  tibble::rownames_to_column(var = "Barcode") %>%
  dplyr::arrange(seurat_clusters)

Cluster[["seurat_clusters"]] <- as.factor(Cluster[["seurat_clusters"]])
write.csv(Cluster, file = file.path(outdir, "seurat_clusters.csv"), quote = FALSE, row.names = FALSE)

cat("All done. Results written to: ", outdir, "\n")
