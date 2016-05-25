sink(file="/tmp/none")
sink("/dev/null")
options(warn=-1)
options(echo=F) 

invisible("EBSeq")
suppressMessages(library("EBSeq"))

args <- commandArgs()
inputfile <- args[6]
outputfile <- args[7]
#PairwisePlots <-args[6]

print(args)

Conds=strsplit(inputfile,split=",")[[1]]


Out=GetPatterns(Conds)


#write.table(Out,file=outputfile,quote=F,col.names=T,row.names=T,sep = "\t")
write.table(Out,file=outputfile,quote=F,col.names=NA,row.names=T,sep = "\t")  ##cms - changed on 11/11/13

