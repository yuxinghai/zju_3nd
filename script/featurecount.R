library("Rsubread")
gtf<- '/data2/zhoulab/yuxinghai/zju/anno/ce11/c_elegans.PRJNA13758.WS258.canonical_geneset.gtf'
name <-c("911","912","913","914","915","916","917","918")
raw_dir<-'/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/02_mapping/sorted'
DESdir <-'/data2/zhoulab/yuxinghai/zju/3nd_analysis/results/03_featurecount'
for (i in name) {
  bsname <-paste(i,"uniq_sort.bam",sep = "_")
  filename <-paste(raw_dir,i,bsname,sep = "/")
  fc_SE <- featureCounts(filename,annot.ext=gtf,isGTFAnnotationFile=T,isPairedEnd=TRUE,strandSpecific
                         =0,nthreads=8)
  write.table(fc_SE$counts,paste(DESdir,paste(i,"featureCounts.txt",sep = "_"),sep = "/"),quote = F,col.names = F,sep = "\t")
  write.table(fc_SE$stat,paste(DESdir,paste(i,"featureStat.log",sep = "_"),sep = "/"),quote = F,col.names = F,sep = "\t")
  
}

