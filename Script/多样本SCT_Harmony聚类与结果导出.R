#20240429更新，添加了单样品的loom文件生成
#s1.2_20240513更新，添加不同聚类cluster导出
#R包加载
library(Seurat)
library(sctransform)
library(ggplot2)
library(harmony)
library(SeuratDisk)
library(stringr)
library(dplyr)
library(tibble)
rm(list=ls())

######参数设置
opt=list()
#res为resolution,值越大聚类个数越多
opt$res=0.3
#筛选基因，表达基因A的最少细胞数量
opt$MinCell=10
#筛选细胞因，细胞中至少表达基因数量
opt$MinFeatures=50
#pca的降维数量
opt$npca=30
#多样本所在的文件夹
opt$indir='path/to/input_dir'
#样本名称
opt$sample=c("10DPA","16DPA","25DPA","35DPA","45DPA")
#结果输出路径
opt$outdir='path/to/output_dir'
#自定义每个cluster的颜色
col = c("0"="#F56867", "1"="#FEB915", "2"="#C798EE","3"="#59BE86","4"="#7495D3","5"="#D1D1D1",
        "6"="#6D1A9C","7"="#15821E","8"="#3A84E6","9"="#70e014","10"="#787878",
        "11"="#DB4C6C","12"="#0430e0","13"="#554236","14"="#AF5F3C","15"="#ff7700","16"="#e00417",
        "17"="#DAB370","18"="#fcfc05","19"="#268785","20"="#09f9f5","21"="#246b93","22"="#cc8e12",
        "23"="#d561dd","24"="#c93f00","25"="#ddd53e","26"="#4aef7b","27"="#e86502","28"="#9ed84e",
        "29"="#39ba30","30"="#6ad157","31"="#8249aa","32"="#99db27","33"="#e07233","34"="#ff523f")


####参数传递
res=opt$res
MinCell=opt$MinCell
MinFeatures=opt$MinFeatures
npca=opt$npca
indir=opt$indir
outdir=paste0(opt$outdir,"/","SCT","_C",MinCell,"_F",MinFeatures,"_P",npca,"_R",res)
if(!dir.exists(outdir)){dir.create(path=outdir,recursive = TRUE)}
dir=c(paste0(indir,"/",opt$sample))
names(dir)=opt$sample
names=opt$sample

####基本分析
#创建seurat对象
expr <- list()
obj <- list()
single<-list()
for(i in 1:length(dir)){
  expr[[i]]=Read10X(dir[i],cell.column = 1)
  obj[[i]]=CreateSeuratObject(counts = expr[[i]],assay = "RNA",min.cells=MinCell,min.features=MinFeatures)
  obj[[i]]@meta.data$orig.ident <-names(dir)[i]
  print(i)
}
object <- merge(x=obj[[1]],y=obj[2:length(obj)])

#数据标准化、高变基因筛选和数据缩放，NFS和SCT是两种方法，选择其一
object <- SCTransform(object,verbose = F,assay = 'RNA')
#PCA降维
object <- RunPCA(object, npcs = 50, verbose = FALSE,reduction.name = "pca") 
#筛选PCA个数，确定npca值，通常选择曲线平滑后的数值
plot_elb <- ElbowPlot(object,ndims = 50)
plot_elb
ggsave(filename = paste(outdir,"/","elbowplot.png",sep=""), plot = plot_elb, width = 10, height = 8)
#不去批次进行聚类
object <- FindNeighbors(object, dims = 1:npca, reduction = "pca")
object <- FindClusters(object, resolution = res, cluster.name = "unintegrated_clusters")
object <- RunUMAP(object, dims = 1:npca, reduction = "pca", reduction.name = "umap.unintegrated")
plot_iu <- DimPlot(object,reduction = "umap.unintegrated",group.by = c("orig.ident","unintegrated_clusters"),pt.size = 1)
plot_iu
ggsave(filename = paste(outdir,"/","umap_initial_unintegrated_clusters.png",sep=""), plot = plot_iu, width = 16, height = 8)
#Harmony2方法进行去批次
object <- RunHarmony(object,reduction="pca", group.by.vars = "orig.ident",reduction.save="integrated.harmony2")
# #CCA方法进行去批次
# options(future.globals.maxSize = 3e+09)
# object <- IntegrateLayers(
#   object = object, method = CCAIntegration,
#   normalization.method = "SCT", new.reduction = "integrated.cca",
#   verbose = FALSE
# )
# #RPCA方法进行去批次
# object <- IntegrateLayers(
#   object = object, method = RPCAIntegration,
#   normalization.method = "SCT", new.reduction = "integrated.rpca",
#   verbose = FALSE
# )
# #harmony方法进行去批次
# object <- IntegrateLayers(
#   object = object, method = HarmonyIntegration,
#   normalization.method = "SCT", new.reduction = "integrated.harmony",
#   verbose = FALSE
# )
#去批次后进行聚类
for (i in c("harmony2")){
  object <- FindNeighbors(object, reduction = paste0("integrated.",i), dims = 1:npca)
  object <- FindClusters(object, resolution = res, cluster.name = paste0(i,"_clusters"))
  object <- RunUMAP(object, reduction = paste0("integrated.",i), dims = 1:npca, reduction.name = paste0("umap.",i))
}
#去批次后作图
for (i in c("harmony2")){
  umapplotorig <- DimPlot(object,reduction = paste0("umap.",i),group.by = "orig.ident",pt.size = 2)
  umapplotorig
  ggsave(filename = paste(outdir,"/","umap_",i,"_orig.ident.png",sep=""), plot = umapplotorig, width = 10, height = 8)
  umapplot <- DimPlot(object,reduction = paste0("umap.",i),group.by = paste0(i,"_clusters"),pt.size = 2,cols = col)
  umapplot
  ggsave(filename = paste(outdir,"/","umap_",i,".png",sep=""), plot = umapplot, width = 10, height = 8)
}

#导出聚类信息
for (j in c("harmony2")){
  cluster_path=paste0(outdir,"/",j)
  if(!dir.exists(cluster_path)){dir.create(path=cluster_path,recursive = TRUE)}
  Cluster=object@meta.data %>%
    select(.,c(paste0(j,"_clusters"))) %>%
    rownames_to_column(.,var = "Barcode") %>%
    arrange(.,paste0(j,"_clusters"))
  Cluster[[paste0(j,"_clusters")]]<- as.factor(Cluster[[paste0(j,"_clusters")]])
  #获得输出文件路径
  for (i in names){
    #获得聚类文件
    filer <- filter(Cluster,grepl(i, Barcode))
    filer$Barcode <- str_replace(filer$Barcode, paste0(i,"_"), "")
    write.csv(filer,file=paste0(cluster_path,"/",i,"_cluster.csv"),quote = F,row.names = F)
  }
}

#发现特征基因
for (i in c("harmony2")){
  Idents(object)=paste0(i,"_clusters")
  #发现特征基因
  object <- PrepSCTFindMarkers(object)
  object.markers <- FindAllMarkers(object, assay = "SCT",only.pos = TRUE)
  write.csv(object.markers,file = paste0(outdir,"/",i,"_marker_gene_onlypos.csv"),row.names = FALSE)
  #拟bulk基因平均表达
  aver.exp <- as.data.frame(AverageExpression(object,assays = "SCT", verbose=F))
  write.csv(aver.exp,file = paste0(outdir,"/",i,"_average_expression.csv"),row.names = TRUE)
}

#保存RDS
saveRDS(object,file = paste0(outdir,"/","object.rds"))
#保存各个样品RDS和loom文件
for (i in names){
  subobject <- subset(object,subset = orig.ident== i)
  colnames(subobject) <- str_replace(colnames(subobject), paste0(i,"_"), "")
  print(paste0("Star save RDS/loom file:",i))
  saveRDS(subobject,file = paste0(outdir,"/",i,"_object.rds"))
  sdata.loom <- as.loom(x = subobject, filename = paste0(outdir,"/",i,"_object.loom"), verbose = FALSE)
  # Always remember to close loom files when done
  sdata.loom$close_all()
}
