sink(file="/tmp/none")
sink("/dev/null")
options(warn=-1)
options(echo=F) 

#invisible("EBSeq")
#suppressMessages(library("EBSeq"))  ##cms - changed on 11/11/13. Don't need EBSeq package.

args <- commandArgs(trailingOnly = T)
inputfile <- args[1]
Idx <- args[2]
outputfile <- args[3]

print(args)

ReadIn=read.table(inputfile,stringsAsFactors=F,header=T, sep="\t")

IndexIn=strsplit(Idx,split=",")[[1]]
Index=as.numeric(IndexIn)

#Mat=data.matrix(ReadIn)
Mat=data.matrix(ReadIn[,-1])  ##cms - changed on 11/11/13. Otherwise pattern IDs are replaced with NA.
rownames(Mat)=ReadIn[,1]  ##cms - also added on 11/11/13 to add pattern IDs back in.

Out=Mat[Index,]


#write.table(Out,file=outputfile,quote=F,col.names=T,row.names=T,sep = "\t")
write.table(Out,file=outputfile,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - changed on 11/11/13

