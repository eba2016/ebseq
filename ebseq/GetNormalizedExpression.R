sink(file="/tmp/none")
sink("/dev/null")
options(warn=-1)
options(echo=F) 

invisible("EBSeq")
suppressMessages(library("EBSeq"))

args <- commandArgs(trailingOnly = T)
inputfile <- args[1]
WhetherSampleName <- args[2]
outputfile <- args[3]
Boxplots<-args[4]
Sizesout <-args[5]

print(args)

if(WhetherSampleName=="y"){
	ReadIn=read.table(inputfile,stringsAsFactors=F,header=T, sep="\t")
	Names=names(ReadIn)[-1]
	}
if(WhetherSampleName=="n"){
	ReadIn=read.table(inputfile,stringsAsFactors=F,header=F, sep="\t")
	Names=paste0("S",1:(length(ReadIn[1,])-1))  ##cms - added 11/11/13
}

GeneMat=do.call(cbind,ReadIn[-1])
rownames(GeneMat)=ReadIn[[1]]
if(WhetherSampleName=="y")colnames(GeneMat)=Names
if(WhetherSampleName=="n")colnames(GeneMat)=Names  ##cms - added 11/11/13

Sizes=MedianNorm(GeneMat)

Data.norm=GetNormalizedMat(GeneMat, Sizes)

#write.table(round(Data.norm,2),file=outputfile,quote=F,col.names=T,row.names=T,sep = "\t")
write.table(round(Data.norm,2),file=outputfile,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - fixed on 11/11/13
pdf(Boxplots)
boxplot(Data.norm,log="y",ylim=c(10^-1,10^6),xlab="Sample",ylab="Normalized Expression")  ##cms - added axis labels on 11/11/13
dev.off()

write.table(Sizes,file=Sizesout,quote=F,col.names=F,row.names=F,sep = "\t")

