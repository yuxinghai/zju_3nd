setwd("/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount")
period <-"h_12"
genename_dir <- paste("gene_name",period,sep="/")
if (!file.exists(genename_dir)){
    dir.create(genename_dir)
}

figure_dir <- paste("figure",period,sep="/")
if (!file.exists(genename_dir)){
    dir.create(figure_dir)
}

DEgene_dir <- paste("DEgene",period,sep="/")
if (!file.exists(DEgene_dir)){
    dir.create(DEgene_dir)
}

SCF_dir <-paste("gene_name","gene_SCF",period,sep="/")
if (!file.exists(SCF_dir)){
    dir.create(SCF_dir)
}


# first 4 sample is control,later 4 is damage
library(DESeq2)
directory <-"/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount"
sampleFiles <- grep("_featureCounts.txt",list.files(directory),value=TRUE)
indx= c(seq(27,32),17,18)
h_12 <- sampleFiles[indx]
#head(sampleFiles)
sampleCondition <-c("control","control","control","control","damage","damage","damage","damage")

sampleTable <- data.frame(sampleName = h_12,
                          fileName = h_12,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c('control','damage'))
dds<-DESeq(ddsHTSeq)
colnames(dds) <-gsub("_featureCounts.txt","",h_12)
# no filter the gene count number 
res<-results(dds)
#write all expression gene name to a file
all_exp <-res[res$baseMean !=0,]
res<-res[order(res$padj),]
mcols(res, use.names=TRUE)
#head(res)

#add gene symbol and filter diff expressed gene
#up:log2FC>0.58, down:log2fc<-0.58
res <-res[complete.cases(res),] 
head(res)
res <-res[ res$padj <0.05,]
res$name <- rownames(res)

new <-as.data.frame(res)
down_gene <- new[new$log2FoldChange<(-0.58),]
#sort by -log2fc
down_gene <-down_gene[with(down_gene,order(-abs(log2FoldChange),padj)),]

up_gene <- new[new$log2FoldChange>log2(1.5),]
up_gene <-up_gene[with(up_gene,order(-abs(log2FoldChange),padj)),]
up_and_down<-new[(new$log2FoldChange<(-0.58)) | (new$log2FoldChange>log2(1.5)), ]
#genename for go term

write.table(up_and_down$name,paste(genename_dir,"sig",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)
write.table(down_gene$name,paste(genename_dir,"down",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)
write.table(up_gene$name,paste(genename_dir,"up",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)


#MAplot
png(paste(figure_dir,"MAplot.png",sep="/"))
plotMA(dds,ylim=c(-12,12),main='wound vs control')
dev.off()


#DISTRIBUTION level of DE 

notsig <-dim(all_exp)[1]-dim(up_and_down)[1]
gene_num <-c(dim(up_gene)[1],notsig,dim(down_gene)[1])
expre_level <-c("Up","Not_change","Down")


plot <-data.frame(gene_num=gene_num,expre_level=expre_level)
library(ggplot2)
ggplot(plot)+geom_bar(aes(x=expre_level,y=gene_num),stat = "identity",
                      position = "dodge",width = 0.5)+labs(y="Gene number")+
  geom_text(aes(x=expre_level,y=gene_num,label=gene_num),vjust = -0.1 )
ggsave(paste(figure_dir,"DE_number.pdf",sep="/"))

#######################ano siggene####################

ano <- read.csv2("/data2/zhoulab/yuxinghai/zju//anno/ce11/ws258mart_export.csv",
                 sep="\t",header =T ,comment.char = "#")
down_gene <-merge(down_gene,ano,by.x="name",by.y="Gene.stable.ID",all.x=T)
#head(down_gene)
up_gene <-merge(up_gene,ano,by.x="name",by.y="Gene.stable.ID",all.x=T)
write.csv(down_gene[,c(1,8,2,3,6,7,9,10,11)],paste(DEgene_dir,"ano_downgene.csv",sep="/"),quote = F,row.names = F)
write.csv(up_gene[,c(1,8,2,3,6,7,9,10,11)],paste(DEgene_dir,"ano_upgene.csv",sep="/"),quote = F,row.names = F)

write.table(down_gene$Gene.name,paste(SCF_dir,"down_Gene_sym",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)
write.table(up_gene$Gene.name,paste(SCF_dir,"up_Gene_sym",sep="/"),sep = "\n",quote = F,col.names = F,row.names = F)
