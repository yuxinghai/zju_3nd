# transient regulated genes
directory <-"/home/yuxh/zheDA/3nd_analysis/03_featurecount"
setwd(directory)
#distance in heatmap(all samples)

# rank sample 
pVal<- 0.01
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
h_2 <- sampleFiles[9:16]
h_2c <-h_2[1:4] # 4
h_2w <-h_2[5:8] # 4
h_8 <- sampleFiles[19:26]
h_8c <-h_8[1:4] # 4
h_8w <-h_8[5:7] # 3
h_24 <- sampleFiles[1:8]
h_24c <-h_24[2:4] # 3
h_24w <-h_24[5:8] # 4
sampleFiles_new <-c(h_2w,h_8w,h_24w,h_2c,h_8c,h_24c)
name <-gsub("_featureCounts.txt","",sampleFiles_new)
# set annotation
wound <-c(rep(2,4),rep(8,3),rep(24,4))
control <-c(rep(2,4),rep(8,4),rep(24,3))
wound_batch <-c(rep("T2",4),rep("T8",3),rep("T24",4))
control_batch <-c(rep("C2",4),rep("C8",4),rep("C24",3))
library(ImpulseDE2)
#distance in heatmap(all samples)
load_data <- function(sampleFiles_new) {
  files <- sampleFiles_new
  tables <- lapply(files, read.table)
  do.call(cbind, tables)
}
pollutantmean <- load_data(sampleFiles_new)
n_list <-seq(1,ncol(pollutantmean))
index= n_list[n_list %% 2 == 0]
inData <-pollutantmean[,index]
colnames(inData) <-name
rownames(inData) <-pollutantmean[,1]
inData <- as.matrix(inData[rowSums(inData)>0,])
head(inData)
design <- data.frame("Sample"=colnames(inData),
                     "Condition"=c(rep("case",length(wound)),rep("control",length(control))),
                     "Time"=c(wound,control),
                     "Batch"=c(wound_batch,control_batch), 
                     row.names=colnames(inData))

impulse_results <- runImpulseDE2(matCountData = inData,
                                 dfAnnotation =design,
                                 boolCaseCtrl = TRUE,
                                 scaNProc = 50,
                                 scaQThres = pVal,
                                 vecConfounders = c("Batch"),
                                 boolIdentifyTransients = TRUE)

write.table(impulse_results$dfImpulseDE2Results[,c(1,3)],file.path(outPath,outFile),row.names=F,col.names=F,quote=F,sep="\t")
library(ComplexHeatmap)
lsHeatmaps <- plotHeatmap(
  objectImpulseDE2       = objectImpulseDE2,
  strCondition           ="case",
  boolIdentifyTransients = TRUE,
  scaQThres              = 0.01)
pdf("impulse2_heatmap.pdf")
draw(lsHeatmaps$complexHeatmapRaw) 
dev.off()

