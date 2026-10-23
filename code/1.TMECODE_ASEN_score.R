## The ASEN score was used to quantify the colocalization of cell types in the same spot. 
## By establishing the ASEN algorithm, we calculate the stable expression for each gene in the spot and estimate the detection probability of a cell cluster based on the number of cell cluster-specific genes expressed. 
## Then, a permutation method generates a background distribution and the detection probability for each cell cluster is calculated, 
## and their product gives the final probability of detecting multiple cell clusters simultaneously within that spot (see Methods)

cal.prob=function(seurat.obj,marker.list,quantile=0.05){
  binary.mat=t(apply(seurat.obj@assays$Spatial@data,1,function(gexpr){
    gexpr1=gexpr[gexpr>0]
    cutoff=quantile(gexpr1,quantile)
    gexpr2=ifelse(gexpr>cutoff,1,0)
  }))
  binary.mat=binary.mat[!is.na(binary.mat[,1]),]
  marker.list=lapply(marker.list,intersect,rownames(binary.mat))
  prob.list=list()
  for(i in 1:length(marker.list))
  {gs=marker.list[[i]]
  total.count=apply(binary.mat[gs,],2,sum)
  bg.dist=t(sapply(1,function(x){
    binary.mat.rand=apply(binary.mat[gs,],1,function(x1){
      sample(x1,length(x1),replace = F)
    })
    y=apply(binary.mat.rand,1,sum)
  }))
  mean.dist=mean(bg.dist)
  sd.dist=sd(bg.dist)
  prob.list[i]=list(pnorm(total.count,mean = mean.dist,sd.dist))
  }
  prob=apply(eval(as.call(c(rbind,prob.list))),2,function(p){prod(p)})
}

TumorST <- readr::read_rds('./BoundaryDefine.rds.gz')
TumorST <- NormalizeData(TumorST)
TME_Signature <- openxlsx::read.xlsx('./TME_DefineTypesDiffMarkers.xlsx')
FourTypesSignature <- lapply(split(TME_Signature,TME_Signature$cluster)[c('SPP1+ Macro','CCL3+ Neutrophil','FAP+ Fibroblasts','Tip Cells')], function(x)x$gene[1:50])

for (x in c('SPP1+ Macro','CCL3+ Neutrophil','FAP+ Fibroblasts','Tip Cells')) {
  TumorST@meta.data[,x] <- apply(TumorST@assays$Spatial@data[rownames(TumorST@assays$Spatial@data) %in% FourTypesSignature[[x]],],2,mean)
}

TumorST@meta.data$'SignatureScore' <- apply(TumorST@meta.data[,c('SPP1+ Macro','Tip Cells','FAP+ Fibroblasts','CCL3+ Neutrophil')],1,mean)
prob <- cal.prob(TumorST,FourTypesSignature,quantile = 0.05)
TumorST$ASEN_Score <- TumorST$SignatureScore*prob
                             
