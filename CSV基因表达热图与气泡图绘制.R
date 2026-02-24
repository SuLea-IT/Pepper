#基础画图
library(readxl)
library(pheatmap)
library(RColorBrewer)
setwd("path/to/input_dir")

# 读取新格式的CSV数据（基因为行，cluster为列）
retu=read.csv("path/to/input_file.csv", header = TRUE, row.names = 1)
print("原始数据:")
print(head(retu))

# 检查数据结构
print("数据维度:")
print(dim(retu))
print("列名:")
print(colnames(retu))

# 数据已经是基因为行，cluster为列的格式
# 直接使用这个数据作为热图数据
heatmap_data <- retu

print("用于热图的数据:")
print(head(heatmap_data))

# 加载ggplot2用于绘制气泡图
library(ggplot2)
library(dplyr)

# 准备气泡图数据 - 将宽格式转换为长格式
library(reshape2)

# 添加基因名称列
bubble_data <- retu
bubble_data$gene <- rownames(bubble_data)

# 将数据从宽格式转换为长格式
bubble_data_long <- melt(bubble_data,
                         id.vars = "gene",
                         variable.name = "cluster",
                         value.name = "expression")

# 重命名cluster列的值，确保是字符格式
bubble_data_long$cluster <- as.character(bubble_data_long$cluster)

# 为气泡图添加必要的列
bubble_data_long <- bubble_data_long %>%
  mutate(
    # 使用expression作为平均表达量和表达百分比的基础
    avg_expression = log2(expression + 1),  # 对表达量取log2
    # 根据表达量计算表达百分比（这里简化处理）
    pct_expressed = pmin(expression / 10, 100),  # 简化的百分比计算
    # cluster已经是正确的格式
    identity = factor(cluster, levels = c("cluster0", "cluster1", "cluster2"))
  )

print("气泡图数据:")
print(head(bubble_data_long))
print("数据维度:")
print(dim(bubble_data_long))


p1 <- pheatmap(heatmap_data,  #要绘制热图的矩阵
               #color = col,
               color = colorRampPalette(c('#ADD8E6', "#FFFFFF", '#FF0000'))(100), #热图色块颜色是从蓝到红分为100个等级
               border = FALSE,
               cellwidth = 30, 
               cellheight = 10,
               border_color = "#FFFFFF",  #热图中每个色块的边框颜色，NA表示无边框
               scale = "row", #“row”按行进行归一化-单个基因不同组织内归一化，"column"表示按列，"none"表示不进行归一化==基因在那就归一化那列或那行
               cluster_rows = FALSE, #是否对行进行聚类
               cluster_cols = TRUE, #是否对列进行聚类
               treeheight_col = 8,
               treeheight_row = 8,
               legend = TRUE, #是否显示图例
               #legend_breaks = c(-0, 4), #设置图例的断点
               #legend_labels = c("0","4"), #设置图例断点处的标签
               show_rownames = TRUE, #是否显示行名
               show_colnames = TRUE, #是否显示列名
               #main ="Female vs Male (Ovary)",
               fontsize_row = 5, #设置标签大小和倾斜角度;#字体大小，可以通过fontsize_row、fontsize_col参数分别设置行列名的字体大小
               fontsize_col = 7,
               angle_col = 0,
               #display_numbers = FALSE,  #是否显示每个色块对应的数值(经归一化后的数值)
               # display_numbers = matrix(ifelse(retu1 > 1, "*", " "), nrow = nrow(retu1)),
               # number_format = "%.2f",  #数值格式，%.2f表示保留小数点后两位;#%.1e表示使用科学计数法并保留小数点后一位
               # number_color = "black",  #设置数值颜色
               # fontsize_number = 4,     #设置数值的字体大小
               #annotation_col = annotation_col, #定义注释栏
               #annotation_colors = ann_colors
) #定义注释栏颜色

pdf("lc-heatmap2.pdf",width = 25,height = 30)
p1
dev.off()
