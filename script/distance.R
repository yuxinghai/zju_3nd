setwd("/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount")
#distance in heatmap(all samples)
deseq2_prepare <-function(name,sampleFiles){
  library(DESeq2)
  sampleCondition <-c("control","control","control","control","wound","wound","wound","wound")

  sampleTable <- data.frame(sampleName = sampleFiles,
                            fileName = sampleFiles,
                            condition = sampleCondition)
  ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                         directory = directory,
                                         design= ~ condition)
  colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','wound'))
  dds<-DESeq(ddsHTSeq)
  name_base <-gsub("_featureCounts.txt","",sampleFiles)
  colnames(dds) <-c(paste0("C",name_base[1:4],sep=""),paste0("W",name_base[5:8],sep=""))
  res<-results(dds)
  sampleCor <- cor(assay(dds) )
  as.matrix( sampleCor )
  sampleCorMatrix <- as.matrix( sampleCor )
  
  return(sampleCorMatrix)
}
library("ComplexHeatmap")
library("RColorBrewer")

distance_plot <- function(sample,method) {
  library("ComplexHeatmap")
  library("RColorBrewer")
  colours = colorRampPalette( brewer.pal(9, "Blues") )(255)
  bs_name <-paste(sample,method,"heatmap_sample_distance.pdf",sep="_")
  pdf(paste("figure",bs_name,sep="/"))
  print(Heatmap(sampleCorMatrix,col =colours ,name = "Correlation",clustering_distance_rows = method,clustering_distance_columns = method,heatmap_legend_param = list(color_bar = "continuous")))
  dev.off()

}
directory <-"/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount"
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)


# 24 h 
library("ComplexHeatmap")
library("RColorBrewer")
h_24 <- sampleFiles[1:8]
colours = colorRampPalette( brewer.pal(9, "Blues") )(255)
sampleCorMatrix <-deseq2_prepare("h_24",h_24)  
distance_plot("h_24","euclidean")
distance_plot("h_24","kendall")
distance_plot("h_24","spearman")
pdf("figure/h_24_k-means2_heatmap_sample_distance.pdf")
Heatmap(sampleCorMatrix,col =colours ,name = "Correlation",km=2,cluster_columns = FALSE,heatmap_legend_param = list(color_bar = "continuous"))
dev.off()

# 2h 
h_2 <- sampleFiles[9:16]
sampleCorMatrix <-deseq2_prepare("h_2",h_2)  
distance_plot("h_2","euclidean")
distance_plot("h_2","kendall")
distance_plot("h_2","spearman")
pdf("figure/h_2_k-means2_heatmap_sample_distance.pdf")
Heatmap(sampleCorMatrix,col =colours ,name = "Correlation",km=2,cluster_columns = FALSE,heatmap_legend_param = list(color_bar = "continuous"))
dev.off()
# 8h
h_8 <- sampleFiles[19:26]
sampleCorMatrix <-deseq2_prepare("h_8",h_8)  
distance_plot("h_8","euclidean")
distance_plot("h_8","kendall")
distance_plot("h_8","spearman")
pdf("figure/h_8_k-means2_heatmap_sample_distance.pdf")
Heatmap(sampleCorMatrix,col =colours ,name = "Correlation",km=2,cluster_columns = FALSE,heatmap_legend_param = list(color_bar = "continuous"))
dev.off()
# 12h
indx= c(seq(27,32),17,18)
h_12 <- sampleFiles[indx]
sampleCorMatrix <-deseq2_prepare("h_12",h_12)  
distance_plot("h_12","euclidean")
distance_plot("h_12","kendall")
distance_plot("h_12","spearman")
pdf("figure/h_12_k-means2_heatmap_sample_distance.pdf")
Heatmap(sampleCorMatrix,col =colours ,name = "Correlation",km=2,cluster_columns = FALSE,heatmap_legend_param = list(color_bar = "continuous"))
dev.off()
