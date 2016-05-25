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
FDR <- args[4]
outputfile <- args[5]
Sort.out<-args[6]
Sort.out.FDR <-args[7]
Sizesout <-args[8]

Conditions=strsplit(CondIn,split=",")[[1]]

if(WhetherSampleName=="y"){
	ReadIn=read.table(inputfile,stringsAsFactors=F,header=T,sep="\t")
	Names=names(ReadIn)[-1]
	}
if(WhetherSampleName=="n"){
	ReadIn=read.table(inputfile,stringsAsFactors=F,header=F,sep="\t")
	Names=paste0("S",1:length(Conditions))
}

if(class(ReadIn[[1]])=="character"){
	GeneMat=do.call(cbind,ReadIn[-1])
	rownames(GeneMat)=ReadIn[[1]]
	colnames(GeneMat)=Names
}
if(class(ReadIn[[1]])=="numeric"){
	GeneMat=data.matrix(ReadIn)
	colnames(GeneMat)=Names
	}


Sizes=MedianNorm(GeneMat)
EBOut=EBTest(Data=GeneMat,Conditions=as.factor(Conditions),sizeFactors=Sizes, maxround=5)
PP=GetPP(EBOut)
PP.sort=sort(PP,decreasing=T)
PP.sort.FDR=PP.sort[which(PP.sort>=1-as.numeric(FDR))]

Data.norm=GetNormalizedMat(GeneMat, Sizes)
FC=PostFC(EBOut)
realFC=FC[[2]]
postFC=FC[[1]]

Mat=cbind(PP, realFC[names(PP)], postFC[names(PP)],Data.norm[names(PP),])
Mat.sort=cbind(PP.sort, realFC[names(PP.sort)], postFC[names(PP.sort)],Data.norm[names(PP.sort),])


if(length(PP.sort.FDR)>1)Mat.sort.FDR=cbind(PP.sort.FDR, realFC[names(PP.sort.FDR)], postFC[names(PP.sort.FDR)],Data.norm[names(PP.sort.FDR),])

if(length(PP.sort.FDR)==1)Mat.sort.FDR=matrix(
		    c(PP.sort.FDR, realFC[names(PP.sort.FDR)], postFC[names(PP.sort.FDR)],Data.norm[names(PP.sort.FDR),])
				    ,nrow=1)

colnames(Mat)=colnames(Mat.sort)=
					  c("PPDE","RealFC","PosteriorFC",colnames(Data.norm))
if(length(PP.sort.FDR)>0)colnames(Mat.sort.FDR)=
			      c("PPDE","RealFC","PosteriorFC",colnames(Data.norm))

##cms - col.names changed below to have a header for the row names. This way the output of RSEM can be joined with EBSeq output
#write.table(round(Mat,2),file=outputfile,quote=F,col.names=T,row.names=T,sep = "\t")
#write.table(round(Mat,2),file=outputfile,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - made correction on 11/11/13
write.table(round(Mat,2),file=outputfile,quote=F,col.names=c("gene_id\tPPDE","RealFC","PosteriorFC",colnames(Data.norm)),row.names=T,sep = "\t")  ##cms - made correction on 11/11/13
#write.table(round(Mat.sort,2),file=Sort.out ,quote=F,col.names=T,row.names=T,sep = "\t")
#write.table(round(Mat.sort,2),file=Sort.out ,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - made correction on 11/11/13
write.table(round(Mat.sort,2),file=Sort.out ,quote=F,col.names=c("gene_id\tPPDE","RealFC","PosteriorFC",colnames(Data.norm)),row.names=T,sep = "\t")  ##cms - made correction on 11/11/13
#if(length(PP.sort.FDR)>0)write.table(round(Mat.sort.FDR,2),file=Sort.out.FDR,quote=F,col.names=T,row.names=T,sep = "\t")
#if(length(PP.sort.FDR)>0)write.table(round(Mat.sort.FDR,2),file=Sort.out.FDR,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - made correction on 11/11/13
if(length(PP.sort.FDR)>0)write.table(round(Mat.sort.FDR,2),file=Sort.out.FDR,quote=F,col.names=c("gene_id\tPPDE","RealFC","PosteriorFC",colnames(Data.norm)),row.names=T,sep = "\t")  ##cms - made correction on 12/17/13; this preserves the header if merged with RSEM output.
write.table(Sizes,file=Sizesout,quote=F,col.names=F,row.names=F,sep = "\t")

