setwd("/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount")
directory <-"/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount"
#distance in heatmap(all samples)
library(DESeq2)
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
#head(sampleFiles)
sampleCondition <-c("control","control","control","control","damage","damage","damage","damage")

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition) 
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','damage'))
dds<-DESeq(ddsHTSeq)
colnames(dds) <-c("911","912","913","914","915","916","917","918")
res<-results(dds)
sampleCor <- cor(assay(dds) )
as.matrix( sampleCor )[ 1:8, 1:8]
sampleCorMatrix <- as.matrix( sampleCor )
#colnames(sampleDistMatrix) <- NULL   
library(ComplexHeatmap)
library("RColorBrewer")
pdf("figure/sample_distance.pdf")
colours = colorRampPalette( brewer.pal(9, "Blues") )(255)
Heatmap(sampleCorMatrix,cluster_rows =F,cluster_columns=F,col =colours ,name = "Correlation",heatmap_legend_param = list(color_bar = "continuous"))
dev.off()
