setwd("/work/result.EN/dr.CCLE/04-mix/09.Expr/")
W5.files = dir("W5/")
dir.create("DAVID")
dir.create("DAVID/W5_Q95")

for(k in 1:length(W5.files)){
	drug = gsub("W5.", "", gsub(".txt", "", W5.files[k]))
	
	mydata = read.table(paste("W5/", W5.files[k], sep=""), header=T, as.is=T)
	colnames(mydata)[1:3] = c("cancer", "gene", "drug")
	mydata = mydata[mydata[,1]!="BRCA" & mydata[,1]!="LGG" & mydata[,1]!="THCA",]
	
	tapply(mydata[,5], mydata[,"tag"], function(u)sum(u)/sqrt(length(u))) -> check
	check = sort(check)
	names(which(check > quantile(check, .95))) -> genes
	gsub(paste("_",drug,sep=""), "", genes) -> genes
	write.table(genes, file=paste("/work/result.EN/dr.CCLE/04-mix/09.Expr/DAVID/W5_Q95/", drug, ".high.Q95.txt", sep=""), row.names=F, quote=F, sep="\t", col.names=F)
	names(which(check < quantile(check, .05))) -> genes
	gsub(paste("_",drug,sep=""), "", genes) -> genes
	write.table(genes, file=paste("/work/result.EN/dr.CCLE/04-mix/09.Expr/DAVID/W5_Q95/", drug, ".low.Q95.txt", sep=""), row.names=F, quote=F, sep="\t", col.names=F)
}

