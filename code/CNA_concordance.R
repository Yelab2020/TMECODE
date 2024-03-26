library(copykat)
library(infercnv)
library(Seurat)
library(dplyr)
MyColorCode <- c('#dc8e97','#e3d1db','#74a893','#ac9141','#5ac6e9','#e5c06e','#7587b1','#c7deef','#e97371','#e1a4c6',
                 '#916ba6','#cb8f82','#7db3af','#d2e0ac','#f3f3b0','#a5d177','#e0bc58','#64abc0','#fab37f','#e98741',
                 '#8fc0dc','#967568','#f2d3ca','#eebd85','#82c785','#edeaa4','#cdaa9f','#ffffff','#bcacd3','#889b5d',
                 '#4e9592','#dbad5f','#64ae79','#64ae79','#ac5092','#f37d74','#4a6eb6','#bbe3ee','#a5d177','#5b9951',
                 '#fbd9ae','#cf622b','#90cdc1','#f4eb9b','#5f9d58','#f5949f','#c1cf9a','#70afab','#a8c375','#8e95c9',
                 '#cf622b','#5384c4','#cacae3','#e45878','#996da8','#3b99d4')

####CNV inferred by CopyKat in scRNAseq
M1026_epi <- readr::read_rds('M1026_epi_SeuObject.rds.gz',compress = 'gz')
exp.rawdata <- GetAssayData(M1026_epi,slot = 'counts')
copykat.test <- copykat(rawmat=exp.rawdata, id.type="S", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=i,
                        distance="euclidean",output.seg="FLASE", plot.genes="TRUE", genome="hg20",n.cores=4)
readr::write_rds(copykat.test,'./copykat.test.rds.gz',compress = 'gz')


####CNV inferred by CNVkit in WES
cnvtable <- read.table('./cnvkit/M1026_sorted_markdup_BQSR.cnr',header = T)
ChrCol <- MyColorCode[1:24] %>% set_names(unique(cnvtable$chromosome))
right_annotation = ComplexHeatmap::rowAnnotation(Chr = cnvtable$chromosome,col = list(Chr = ChrCol))

pdf('./M1026_WES_cnv.pdf',width = 3,height = 10)
lapply(unique(cnvtable$chromosome), function(x){
  tmp <- cnvtable[cnvtable$chromosome %in% x,]
  right_annotation = ComplexHeatmap::rowAnnotation(Chr = tmp$chromosome,col = list(Chr = ChrCol[x]))
  
  cnvmat <- as.matrix(tmp$log2) %>% set_rownames(tmp$chromosome)
  ComplexHeatmap::Heatmap(cnvmat,cluster_columns = F,cluster_rows = F,
                          show_row_names = F,show_heatmap_legend = FALSE,
                          col = colorRamp2(c(-2, 0, 2), c( "#73b0d2", "white", "#ef926d")))
})
dev.off()


####CNV inferred by inferCNV in inferCNV
cnv_table <- read.table('./InferCNV/M1026/infercnv.observations.txt')

gene_order_file <- read.table('./gencode_v38_gene_pos.txt') %>% column_to_rownames('V1')

Chr <- data.frame(Gene = rownames(cnv_table),Chr = gene_order_file[rownames(cnv_table),]$V2)


pdf('./M1026_inferCNV_cnv.pdf',height = 5,width = 5)
lapply(unique(Chr$Chr), function(x){
  ComplexHeatmap::Heatmap(as.matrix(cnv_table)[Chr[Chr$Chr %in% x,]$Gene,],show_heatmap_legend = F,
                          cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F)
  
})
dev.off()
