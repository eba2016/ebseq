sink(file="/tmp/none")
sink("/dev/null")
options(warn=-1)
options(echo=F) 

invisible("EBSeq")
suppressMessages(library("EBSeq"))

args <- commandArgs(trailingOnly = T)
inputfile <- args[1]
outputfile <- args[2]

print(args)

a1=read.csv(inputfile,stringsAsFactors=F,header=F, sep="\t")
Ng=GetNg(a1[[1]],a1[[2]])
Ig=Ng$IsoformNgTrun




write.table(Ig,file=outputfile,quote=F,col.names=F,row.names=F,sep = "\t")

