library(ggplot2)

# ========== 函数：绘制某基因在两个方向的拟时序曲线 ========== #
plot_gene_two_directions <- function(cds, gene_id, gene_sym,
                                     states = c(2, 3),
                                     min_expr = 0.0,
                                     cluster_cols = NULL,
                                     point_size = 2,
                                     span = 0.6) {
  pd <- pData(cds)
  stopifnot(all(c("Pseudotime","State","cluster_stage3") %in% colnames(pd)))
  
  expr_v <- as.numeric(exprs(cds[gene_id, ]))
  df <- data.frame(
    Pseudotime = pd$Pseudotime,
    Expr       = expr_v,
    State      = pd$State,
    Cluster    = pd$cluster_stage3,
    stringsAsFactors = FALSE
  )
  
  # 过滤：只保留目标 State、有效时间点、表达超过阈值
  df <- df[df$State %in% states & !is.na(df$Pseudotime) & df$Expr >= min_expr, , drop = FALSE]
  
  # 定义方向（用于线型）
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

# ========== 批量绘制并保存 ========== #
for (gene_id in names(genes_to_plot)) {
  gene_sym <- genes_to_plot[gene_id]
  min_e    <- min_expr_tbl[gene_sym]
  
  p_gene <- plot_gene_two_directions(
    cds          = mycds_var,       # 你的 CellDataSet
    gene_id      = gene_id,
    gene_sym     = gene_sym,
    states       = c(2, 3),           # 换成你实际想对比的两个 State
    min_expr     = min_e,
    cluster_cols = base_col,
    point_size   = 2,
    span         = 0.6
  )
  
  ggsave(filename = file.path(opt$outdir, paste0(gene_sym, "_branched.pdf")),
         plot = p_gene, width = 6, height = 4.5, device = cairo_pdf)
}
