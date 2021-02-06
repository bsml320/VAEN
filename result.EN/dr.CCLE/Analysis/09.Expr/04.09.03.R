setwd("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/04-mix/09.Expr")
dir.create("W5")
files = dir()
files = files[grep("W5", files)]
files = files[grep("txt", files)]

expr.mat = read.delim(files[1], as.is=T)
drugs = unique(expr.mat[,3])
drugs[which(drugs=="X17.AAG")] = "17-AAG"
gsub("\\.", "-", drugs) -> new.drugs
drugs = new.drugs

#################################################################################################

for(k in 1:length(files)){
	expr.mat = read.delim(files[k], as.is=T)
	expr.mat[which(expr.mat[,3]=="X17.AAG"),3] = "17-AAG"
	ss = expr.mat[,3]
	expr.mat[,3] = gsub("\\.", "-", ss)
	expr.mat$adjp = p.adjust(expr.mat[,4], method="BH")
	grep("ENSG", expr.mat[,1]) -> ii
	if(length(ii) > 0)expr.mat = expr.mat[-ii, ]
	dat.mat = expr.mat[expr.mat[,4] < 0.05, ]
	dat.mat = expr.mat
	apply(dat.mat, 1, function(u)paste(u[2], u[3], sep="_")) -> tag
	dat.mat = cbind(dat.mat, tag)
	
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		skcm.expr.mat = dat.mat[dat.mat[,3]==drug, ]
		if(k==1){
			write.table(skcm.expr.mat, file=paste("W5/",drugs[kdrug], ".txt", sep=""), row.names=F, col.names=T, sep="\t", quote=F)
		} else {
			write.table(skcm.expr.mat, file=paste("W5/",drugs[kdrug], ".txt", sep=""), row.names=F, col.names=F, sep="\t", quote=F, append=T)
		}
	}
	cat(files[k], "\n", sep="")
}

#################################################################################################
