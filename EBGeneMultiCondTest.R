sink(file="/tmp/none")
sink("/dev/null")
options(warn=-1)
options(echo=F) 

invisible("EBSeq")
suppressMessages(library("EBSeq"))

args <- commandArgs(trailingOnly = T)
inputfile <- args[1]
WhetherSampleName <- args[2]
CondIn <- args[3]
PatternFile <- args[4]
outputfile <- args[5]
MAP.out<-args[6]
Sizesout <-args[7]

#write.table(PatternFile,outputfile,quote=F,col.names=T,row.names=T,sep = "\t")

#write.table(args,outputfile,quote=F,col.names=T,row.names=T,sep = "\t")

Conditions=strsplit(CondIn,split=",")[[1]]

if(WhetherSampleName=="y"){
	ReadIn=read.table(inputfile,stringsAsFactors=F,header=T, sep="\t")
	Names=names(ReadIn)[-1]
	}
if(WhetherSampleName=="n"){
	ReadIn=read.table(inputfile,stringsAsFactors=F,header=F, sep="\t")
	Names=paste0("S",1:length(Conditions))
}

#PatternIn=read.table(PatternFile,stringsAsFactors=F,header=T,sep="\t")
PatternIn=read.table(PatternFile,stringsAsFactors=F,header=T,sep="\t",row.names=1) ##cms - changed on 11/15/13

if(class(ReadIn[[1]])=="character"){
	GeneMat=do.call(cbind,ReadIn[-1])
	rownames(GeneMat)=ReadIn[[1]]
	colnames(GeneMat)=Names
}
if(class(ReadIn[[1]])=="numeric"){
	GeneMat=data.matrix(ReadIn)
	colnames(GeneMat)=Names
	}

Patterns=data.matrix(PatternIn)


Sizes=MedianNorm(GeneMat)
EBOut=EBMultiTest(Data=GeneMat,Conditions=as.factor(Conditions),
			AllParti=Patterns,sizeFactors=Sizes, maxround=5)
PPout=GetMultiPP(EBOut)
MultiPP=PPout$PP
MultiMAP=PPout$MAP
Data.norm=round(GetNormalizedMat(GeneMat, Sizes),2)

Mat=cbind(MultiMAP,Data.norm[names(MultiMAP),])

colnames(Mat)=
c("MAP",Names)
options(warn=-1)

#write.table(round(MultiPP,2),file=outputfile,quote=F,col.names=T,row.names=T,sep = "\t")
write.table(round(MultiPP,2),file=outputfile,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - changed on 11/15/13
#write.table(Mat,file=MAP.out ,quote=F,col.names=T,row.names=T,sep = "\t")
write.table(Mat,file=MAP.out ,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - changed on 11/15/13
write.table(Sizes,file=Sizesout,quote=F,col.names=F,row.names=F,sep = "\t")

