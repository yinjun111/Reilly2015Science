args<-commandArgs(TRUE)
infile<-args[1]
outfile<-args[2]

data<-read.table(infile,header=T,row.names=1,sep="\t",check.names =F,quote="\"")
data.pvalue<-apply(data,1,function(x) {binom.test(x[2],x[1],x[3]/x[1],,alternative="greater")$p.value})
data.qvalue<-p.adjust(data.pvalue,method="BH")

#data.fc<-log2(data[,1]/data[,2]/(data[,3]/data[,4]))

data.pvalue.pm<-data[,5]
data.qvalue.pm<-p.adjust(data.pvalue.pm,method="BH")

result<-cbind(data,data.qvalue.pm,data.pvalue,data.qvalue)
colnames(result)<-c(colnames(data),"BH Perm P","Binom P","BH Binom P")

write.table(file=outfile,result,sep="\t",quote=F,col.names=NA)
