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

#### 1. Mfuzz time-series clustering ####
#1).create object：
eset <- new("ExpressionSet",exprs = data.matrix(DataMean.matrix))
#eset <- filter.std(eset,min.std=0)
#2).normalization
eset <- Mfuzz::standardise(eset)
#==omit NA
eset <- Mfuzz::filter.NA(eset, thres=0.25)
#3).clustering
#  estimate best m value
m <- Mfuzz::mestimate(eset)

lapply(2:10, function(i){
  # clusteing
  cl <- mfuzz(eset, c = i, m = m)
  cl$size
  #4). visualization
  colname <- c('Normal','Polyps','I','II','III','IV')
  set.seed(123)
  pdf(file = paste0("1.1_Cluster_C",i,".pdf"),height = 20,width =18)
  mfuzz.plot(eset,cl,mfrow=c(4,3),
             time.labels = colname,
             new.window= FALSE)+
    theme(axis.text.x = element_text(angle = 60))
  
  TimeCluster_sub <- cl$cluster %>% as.data.frame()
  colnames(TimeCluster_sub)[1] <- "cluster"
  # check membership
  TimeCluster_sub$membership <- cl$membership
  TimeCluster_sub$gene <- rownames(sub)
  #Trmt61a,Myc,Trmt6：C2
  write.csv(TimeCluster_sub,paste0("2.Cluster_sub",i,".csv"))
  write.table(cl$size,paste0("3.Cluster_sub",i,"_size.txt"))
  dev.off()
})
