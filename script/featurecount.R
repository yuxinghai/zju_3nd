library("Rsubread")
gtf<- '/data2/zhoulab/yuxinghai/zju/anno/ce11/c_elegans.PRJNA13758.WS258.canonical_geneset.gtf'
name <-c("903","905","907","909","911","913","915","917","923","935","937","939","941","943","945","947","904","906","908","910","912","914","916","918","924","936","938","940","942","944","946","948")
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

