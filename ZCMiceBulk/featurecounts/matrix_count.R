library(reshape2)

# generate count matrix
a=read.csv('featurecounts/gene.count',header=F,sep="\t")

colnames(a)=c('sample','gene','counts')

counts=dcast(a,formula=gene~sample)

write.table(counts,file="featurecounts/gene_count.csv",sep="\t",quote=FALSE,row.names=FALSE)

# generate tpm matrix
a=read.csv('featurecounts/gene.tpm',header=F,sep="\t")

colnames(a)=c('sample','gene','tpm')

counts=dcast(a,formula=gene~sample)

write.table(counts,file="featurecounts/gene_tpm.csv",sep="\t",quote=FALSE,row.names=FALSE)


# 链接：https://www.jianshu.com/p/3352bfd404f3
