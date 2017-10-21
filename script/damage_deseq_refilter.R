setwd("/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount")
# first 3 sample is control,later 3 is damage
library(DESeq2)
directory <-"/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount"
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

write.table(up_and_down$name,"gene_name/sig",sep = "\n",quote = F,col.names = F,row.names = F)
write.table(down_gene$name,"gene_name/down",sep = "\n",quote = F,col.names = F,row.names = F)
write.table(up_gene$name,"gene_name/up",sep = "\n",quote = F,col.names = F,row.names = F)

write.csv(down_gene[,c(7,1,2,5,6)],"DEgene/damage_vs_control_downgene.csv",quote = F,row.names = F)
write.csv(up_gene[,c(7,1,2,5,6)],"DEgene/damage_vs_control_upgene.csv",quote = F,row.names = F)

#MAplot
png("figure/damage_vs_control_MAplot.png")
plotMA(dds,ylim=c(-12,12),main='damage vs control')
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
ggsave("figure/Gene_number_in_DE_level.pdf")

#######################ano siggene####################

ano <- read.csv2("/data2/zhoulab/yuxinghai/zju//anno/ce11/ws258mart_export.csv",
                 sep="\t",header =T ,comment.char = "#")
down_gene <-new <-merge(down_gene,ano,by.x="name",by.y="Gene.stable.ID")
#head(down_gene)
up_gene <-new <-merge(up_gene,ano,by.x="name",by.y="Gene.stable.ID")
write.csv(down_gene[,c(1,8,2,3,6,7,9,10,11)],"DEgene/ano_damage_vs_control_downgene.csv",quote = F,row.names = F)
write.csv(up_gene[,c(1,8,2,3,6,7,9,10,11)],"DEgene/ano_damage_vs_control_upgene.csv",quote = F,row.names = F)

write.table(down_gene$Gene.name,"gene_name/gene_SCF/down_Gene_sym",sep = "\n",quote = F,col.names = F,row.names = F)
write.table(up_gene$Gene.name,"gene_name/gene_SCF/up_Gene_sym",sep = "\n",quote = F,col.names = F,row.names = F)
