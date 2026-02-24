# ===================== 辣椒单细胞 × 空间转录组一体化流程（含UMAP、无回投版 + 自定义调色板） =====================
# v2025-08-23.12
# - 递归从 base_dir 读取 10x (.gz) 三件套并合并
# - 阶段一/二/三聚类；输出各阶段 RDS
# - 输出 UMAP（stage1/2/3），使用自定义调色板 col_primary，但不做任何空间回投
# ===============================================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(dplyr)
  library(Matrix)
  library(ggplot2)   # for ggsave + DimPlot
})

options(future.globals.maxSize = 9e12)
set.seed(1234)

# ===================== 自定义调色板 =====================
col_primary <- c(
  "0"="#F56867", "1"="#FEB915", "2"="#C798EE", "3"="#59BE86", "4"="#7495D3", "5"="#D1D1D1",
  "6"="#6D1A9C", "7"="#15821E", "8"="#3A84E6", "9"="#70E014", "10"="#787878", "11"="#DB4C6C",
  "12"="#0430E0", "13"="#554236", "14"="#AF5F3C", "15"="#FF7700", "16"="#E00417", "17"="#DAB370",
  "18"="#FCFC05", "19"="#268785", "20"="#09F9F5", "21"="#246B93", "22"="#CC8E12", "23"="#D561DD",
  "24"="#C93F00", "25"="#DDD53E", "26"="#4AEF7B", "27"="#E86502", "28"="#9ED84E", "29"="#39BA30",
  "30"="#6AD157", "31"="#8249AA", "32"="#99DB27", "33"="#E07233", "34"="#FF523F", "35"="#D80FC1",
  "36"="#931635", "37"="#280F7A", "38"="#9E7A7A", "39"="#1A918F", "40"="#FF69B4", "41"="#8A2BE2",
  "42"="#00FFFF", "43"="#7FFF00", "44"="#FFD700", "45"="#FF4500", "46"="#DA70D6", "47"="#00CED1",
  "48"="#FF6347", "49"="#4682B4", "50"="#DC143C", "51"="#00FF7F", "52"="#B22222", "53"="#FF8C00",
  "54"="#9932CC", "55"="#8B0000", "56"="#E9967A", "57"="#48D1CC", "58"="#C71585", "59"="#191970",
  "60"="#00BFFF", "61"="#5F9EA0", "62"="#D2691E", "63"="#FF7F50", "64"="#6495ED", "65"="#FFF8DC",
  "66"="#DC82DC", "67"="#2F4F4F", "68"="#4B0082", "69"="#FFE4B5", "70"="#A0522D", "71"="#7B68EE",
  "72"="#556B2F", "73"="#CD853F", "74"="#FFDEAD", "75"="#800000", "76"="#BA55D3", "77"="#32CD32",
  "78"="#3CB371", "79"="#6B8E23", "80"="#8B4513", "81"="#708090", "82"="#ADFF2F", "83"="#8FBC8F",
  "84"="#BDB76B", "85"="#696969", "86"="#FFB6C1", "87"="#FFA500", "88"="#40E0D0", "89"="#1E90FF",
  "90"="#F0E68C", "91"="#B0C4DE", "92"="#9370DB", "93"="#20B2AA", "94"="#8A0303", "95"="#6A5ACD",
  "96"="#9ACD32", "97"="#4B5320", "98"="#FF1493", "99"="#8FBC8B", "100"="#483D8B", "101"="#2E8B57",
  "102"="#B8860B", "103"="#A9A9A9", "104"="#006400", "105"="#8B008B", "106"="#9400D3", "107"="#228B22",
  "108"="#FF00FF", "109"="#DAA520", "110"="#7CFC00", "111"="#ADD8E6", "112"="#F08080", "113"="#E0FFFF",
  "114"="#FAFAD2", "115"="#90EE90", "116"="#FFA07A", "117"="#87CEFA", "118"="#778899", "119"="#FFFFE0",
  "120"="#00FF00", "121"="#66CDAA", "122"="#00FA9A", "123"="#FFE4E1", "124"="#FFEFD5", "125"="#FFDAB9",
  "126"="#FFC0CB", "127"="#DDA0DD", "128"="#B0E0E6", "129"="#800080", "130"="#FF0000", "131"="#BC8F8F",
  "132"="#4169E1", "133"="#FA8072", "134"="#F4A460", "135"="#C0C0C0", "136"="#87CEEB"
)

# 根据分组标签返回命名调色板：优先按名字精确匹配；否则顺序分配；不够则循环
.palette_for_groups <- function(labels_chr) {
  labs <- sort(unique(labels_chr))
  if (length(names(col_primary)) > 0 && all(labs %in% names(col_primary))) {
    pal <- col_primary[labs]
  } else {
    pool <- unname(col_primary)
    if (length(pool) < length(labs)) {
      message(sprintf("注意：颜色数量不足（需要%d，提供%d），将循环使用。", length(labs), length(pool)))
      pool <- rep(pool, length.out = length(labs))
    }
    pal <- pool[seq_along(labs)]
    names(pal) <- labs
  }
  pal
}

# ===================== 配置 =====================
config <- list(
  base_dir    = "path/to/input_dir",            # 放样本文件夹的根目录
  outdir_base = "path/to/output_dir", # 输出目录
  stage1 = list(
    method = "SCT",        # "SCT" 或 "LogNorm"
    resolution = 0.3,
    min_features = 50,
    npca = 30
  ),
  stage2 = list(
    refine_plan = c("0"=0.3, "1"=0.3)   # 为空则跳过
  ),
  stage3 = list(
    reencode = TRUE
  )
)

# ===================== 工具函数 =====================
.ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE)

# 仅识别 .gz 三件套（与你的数据一致）
.is_10x_dir <- function(d) {
  has_mtx  <- file.exists(file.path(d, "matrix.mtx.gz"))
  has_bar  <- file.exists(file.path(d, "barcodes.tsv.gz"))
  has_feat <- file.exists(file.path(d, "features.tsv.gz"))
  has_mtx && has_bar && has_feat
}

# 递归寻找所有 10x 目录
.find_10x_dirs <- function(base_dir) {
  cand <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  dirs <- cand[ vapply(cand, .is_10x_dir, logical(1)) ]
  unique(dirs)
}

# 读取 10x 目录为 Seurat 对象
.read_10x_to_seurat <- function(d, sample_name, min_features = 50) {
  m <- tryCatch(Read10X(data.dir = d), error = function(e) NULL)
  if (is.null(m)) {
    warning("读取失败，跳过：", d); return(NULL)
  }
  if (is.list(m)) {
    if ("Gene Expression" %in% names(m)) m <- m[["Gene Expression"]] else m <- m[[1]]
  }
  obj <- CreateSeuratObject(counts = m, project = sample_name, min.features = min_features)
  obj <- RenameCells(obj, add.cell.id = sample_name)  # 确保跨样本唯一
  obj$orig.ident <- sample_name
  DefaultAssay(obj) <- "RNA"
  obj
}

# 根据已有 PCA 维度返回 1:k
.dims <- function(obj, npca) {
  emb <- tryCatch(Embeddings(obj, "pca"), error=function(e) NULL)
  k <- if (is.null(emb)) npca else min(npca, ncol(emb))
  seq_len(max(1, k))
}

# 主流程（不做空间回投）
run_pipeline <- function(obj, method="SCT", npca=30, resolution=0.3, diag_dir=".") {
  DefaultAssay(obj) <- "RNA"
  if (method == "SCT") {
    obj <- SCTransform(obj, verbose = FALSE)
  } else {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
  }
  obj <- RunPCA(obj, npcs = npca, verbose = FALSE)
  
  # 剔除 PCA 含 NA 的细胞
  pc <- Embeddings(obj, "pca")
  bad <- rownames(pc)[rowSums(is.na(pc)) > 0]
  if (length(bad) > 0) {
    .ensure_dir(diag_dir)
    obj <- subset(obj, cells = setdiff(colnames(obj), bad))
    obj <- RunPCA(obj, npcs = npca, verbose = FALSE)
  }
  
  obj <- FindNeighbors(obj, dims = .dims(obj, npca))
  # 清理旧聚类列
  old <- grep("snn_res|seurat_clusters", colnames(obj@meta.data), value=TRUE)
  if (length(old)>0) for (cc in old) obj[[cc]] <- NULL
  obj <- FindClusters(obj, resolution = resolution)
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) obj$seurat_clusters <- Idents(obj)
  Idents(obj) <- "seurat_clusters"
  obj
}

# —— 生成（或复用）UMAP并出图（使用自定义调色板）；返回更新后的对象 ——
.ensure_umap_and_save <- function(obj, group.by, out_png, npca) {
  if (!("umap" %in% names(obj@reductions))) {
    set.seed(1234)
    obj <- RunUMAP(obj, dims = .dims(obj, npca), verbose = FALSE)
  }
  # 生成按标签命名的调色板
  labs_chr <- as.character(obj@meta.data[, group.by, drop = TRUE])
  pal <- .palette_for_groups(labs_chr)
  p <- DimPlot(obj, reduction = "umap", group.by = group.by, label = TRUE) +
    scale_color_manual(values = pal, drop = FALSE) +
    NoLegend()
  ggsave(out_png, p, width = 8, height = 6, dpi = 300)
  message("UMAP 已保存：", out_png)
  return(obj)
}

# ===================== 主流程 =====================
.out <- config$outdir_base; .ensure_dir(.out)
diag_dir <- file.path(.out, "diagnostics"); .ensure_dir(diag_dir)

# 1) 扫描并读取所有 10x 矩阵目录
base_dir <- config$base_dir
if (!dir.exists(base_dir)) stop("base_dir 不存在：", base_dir)

tenx_dirs <- .find_10x_dirs(base_dir)
if (length(tenx_dirs) == 0) stop("在 base_dir 下未发现任何 10x 矩阵目录：", base_dir)

# 以“目录名”为样本名
seurat_list <- list()
for (d in tenx_dirs) {
  d_norm <- normalizePath(d, winslash = "/", mustWork = FALSE)
  sample_name <- basename(d_norm)
  
  message(sprintf("读取 10x 目录：%s  -> sample=%s", d, sample_name))
  obj <- .read_10x_to_seurat(d, sample_name, min_features = config$stage1$min_features)
  if (!is.null(obj) && ncol(obj) > 0) {
    # 防重名
    if (sample_name %in% names(seurat_list)) {
      k <- sum(startsWith(names(seurat_list), paste0(sample_name, "_"))) + 1
      sample_name <- paste0(sample_name, "_", k)
      obj <- RenameCells(obj, add.cell.id = sample_name)
      obj$orig.ident <- sample_name
    }
    seurat_list[[sample_name]] <- obj
  }
}
if (length(seurat_list) == 0) stop("未成功读取任何样本。")

# 2) 合并
merged <- Reduce(function(a,b) merge(a, b), seurat_list)
merged@graphs <- list()
message("已合并样本：", paste(unique(merged$orig.ident), collapse = ", "))

# 3) 阶段一：主聚类
merged <- run_pipeline(
  merged,
  method     = config$stage1$method,
  npca       = config$stage1$npca,
  resolution = config$stage1$resolution,
  diag_dir   = diag_dir
)
saveRDS(merged, file.path(config$outdir_base, "stage1_merged_clustered.rds"))

# UMAP（stage1）
merged <- .ensure_umap_and_save(
  merged, "seurat_clusters",
  file.path(config$outdir_base, "stage1_umap.png"),
  npca = config$stage1$npca
)

# 4) 阶段二：子聚类（可选）
if (!is.null(config$stage2$refine_plan) && length(config$stage2$refine_plan) > 0) {
  merged$cluster_refined <- as.character(merged$seurat_clusters)
  for (parent_cluster in names(config$stage2$refine_plan)) {
    if (!parent_cluster %in% levels(merged$seurat_clusters)) {
      warning("refine_plan 中的父簇不存在：", parent_cluster); next
    }
    Idents(merged) <- "seurat_clusters"
    sub_obj <- subset(merged, idents = parent_cluster)
    if (ncol(sub_obj) < 10) { warning("子集细胞数过少，跳过父簇 ", parent_cluster); next }
    
    sub_obj <- run_pipeline(
      sub_obj,
      method     = config$stage1$method,
      npca       = config$stage1$npca,
      resolution = config$stage2$refine_plan[[parent_cluster]],
      diag_dir   = diag_dir
    )
    
    sub_lab <- paste0(parent_cluster, "_", as.character(sub_obj$seurat_clusters))
    merged$cluster_refined[Cells(sub_obj)] <- sub_lab
  }
  merged$cluster_refined <- factor(merged$cluster_refined)
  Idents(merged) <- "cluster_refined"
  saveRDS(merged, file.path(config$outdir_base, "stage2_subclustered.rds"))
  
  # UMAP（stage2：复用同一 UMAP 嵌入，只换分组列）
  merged <- .ensure_umap_and_save(
    merged, "cluster_refined",
    file.path(config$outdir_base, "stage2_umap.png"),
    npca = config$stage1$npca
  )
}

# 5) 阶段三：重编码（0..N-1）
if (isTRUE(config$stage3$reencode)) {
  label_col <- if ("cluster_refined" %in% colnames(merged@meta.data)) "cluster_refined" else "seurat_clusters"
  merged$cluster_reencoded <- as.integer(as.factor(merged[[label_col]][,1])) - 1L
  Idents(merged) <- "cluster_reencoded"
  saveRDS(merged, file.path(config$outdir_base, "stage3_reencoded.rds"))
  
  # UMAP（stage3）
  merged <- .ensure_umap_and_save(
    merged, "cluster_reencoded",
    file.path(config$outdir_base, "stage3_umap.png"),
    npca = config$stage1$npca
  )
}
message("流程执行完毕：已输出 UMAP（使用自定义调色板，不做空间回投）。")
