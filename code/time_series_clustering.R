library(Mfuzz)
library(magrittr)
library(clusterProfiler)
library(pheatmap)
library(ggplot2)

##### TimeCluster #####
##Input file: 
#PerList: cell abundance matrix, cell type as rows and sample as columns
#StageList: stage information of samples
PerList <- readr::read_rds('./PerList.rds.gz')
StageList <- readr::read_rds('./StageList.rds.gz')

DataMean.matrix <- lapply(c('Normal','Polyps','I','II','III','IV'), function(x){
  Sample <- StageList[StageList$Stage %in% x,]$Sample
  tmp <- PerList[,as.character(Sample)]
  Mean <- apply(tmp,1,mean) %>% as.data.frame() %>% set_colnames(x)
  return(Mean)
}) %>% bind_cols()

#### 1. Mfuzz 时间序列表达模式聚类分析 ####
#=输入为 FPKM
#1).建立对象：:write_rds(DataMean.matrix,"/work/sjt/Project/Collaboration/LHB/Trmt61/tRNA_AlkB_0,6,18,48h/ProcessedData/FPKM_CodonMean.rds.gz",compress = "gz")

eset <- new("ExpressionSet",exprs = data.matrix(DataMean.matrix))
# 根据标准差去除样本间差异太小的基因
#eset <- filter.std(eset,min.std=0)
#2).标准化
eset <- Mfuzz::standardise(eset)
#==去掉 NA值
eset <- Mfuzz::filter.NA(eset, thres=0.25)
#3).聚类
# 聚类个数

#  评估出最佳的m值
m <- Mfuzz::mestimate(eset)

lapply(2:10, function(i){
  # 聚类
  cl <- mfuzz(eset, c = i, m = m)
  # 查看每个cluster中的基因个数
  cl$size
  #4). 可视化
  colname <- c('Normal','Polyps','I','II','III','IV')
  set.seed(123)
  pdf(file = paste0("1.1_Cluster_C",i,".pdf"),height = 20,width =18)
  mfuzz.plot(eset,cl,mfrow=c(4,3),#2行3列
             time.labels = colname,
             new.window= FALSE)+
    theme(axis.text.x = element_text(angle = 60))
  
  TimeCluster_sub <- cl$cluster %>% as.data.frame()
  colnames(TimeCluster_sub)[1] <- "cluster"
  # 查看基因和cluster之间的membership
  TimeCluster_sub$membership <- cl$membership
  TimeCluster_sub$gene <- rownames(sub)
  #Trmt61a,Myc,Trmt6：C2
  write.csv(TimeCluster_sub,paste0("2.Cluster_sub",i,".csv"))
  write.table(cl$size,paste0("3.Cluster_sub",i,"_size.txt"))
  dev.off()
})
