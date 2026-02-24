# ===================== 导出每个cluster的marker为CSV（含合并版） =====================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

# -------- 配置（按需修改） --------
INPUT_RDS <- "path/to/input_file.rds"   # 你的Seurat对象 .rds
OUT_DIR   <- dirname(INPUT_RDS)                 # 输出目录（默认与输入同目录）
ONLY_POS  <- TRUE                               # 只要上调marker
MIN_PCT   <- 0.10                               # FindAllMarkers: min.pct
LOGFC_TH  <- 0.25                               # FindAllMarkers: logfc.threshold
TEST_USE  <- "wilcox"                           # 统计方法

# -------- 读取对象 --------
obj <- readRDS(INPUT_RDS)

# 若未设置Idents，尝试用 seurat_clusters
if (is.null(Idents(obj)) && !is.null(obj$seurat_clusters)) {
  Idents(obj) <- obj$seurat_clusters
}

# -------- 计算marker --------
markers <- FindAllMarkers(
  obj,
  only.pos        = ONLY_POS,
  min.pct         = MIN_PCT,
  logfc.threshold = LOGFC_TH,
  test.use        = TEST_USE
)

# -------- 列名标准化到所需格式 --------
# Seurat版本差异：可能是 avg_logFC 或 avg_log2FC
if ("avg_logFC" %in% names(markers) && !"avg_log2FC" %in% names(markers)) {
  markers <- markers %>% rename(avg_log2FC = avg_logFC)
}
# 若没有 gene 列，则用行名补齐
if (!"gene" %in% names(markers)) {
  markers$gene <- rownames(markers)
}

# 只保留并按顺序排列所需列
need_cols <- c("p_val","avg_log2FC","pct.1","pct.2","p_val_adj","cluster","gene")
markers <- markers %>% select(all_of(need_cols))

# 合并版：在每个cluster内按 avg_log2FC 降序，并给出分簇内ID
markers_merged <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(ID = row_number()) %>%
  ungroup() %>%
  select(ID, everything())

# 写出合并版
merged_path <- file.path(OUT_DIR, "markers_all_clusters.csv")
write_csv(markers_merged, merged_path)
message("已写出合并版: ", merged_path)
