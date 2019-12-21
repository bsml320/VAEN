setwd("/work/data/ssGSEA")

cancer.types = dir("/work/data/TCGA/")

#> cancer.types
# [1] "ACC"  "BLCA" "BRCA" "CESC" "CHOL" "COAD" "DLBC" "ESCA" "GBM"  "HNSC"
#[11] "KICH" "KIRC" "KIRP" "LAML" "LGG"  "LIHC" "LUAD" "LUSC" "MESO" "OV"
#[21] "PAAD" "PCPG" "PRAD" "READ" "SARC" "SKCM" "STAD" "TGCT" "THCA" "THYM"
#[31] "UCEC" "UCS"  "UVM"


library(GSVA)
load("/work/data/ssGSEA/c2.all.v5.1.symbols.RData")

for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	original.log2.RPKM = read.delim(paste("/work/data/TCGA/",cancer,"/HiSeqV2", sep=""), as.is=T)
	apply(original.log2.RPKM[,-1],1,sum) -> rowCheck
	non0.log2.RPKM = original.log2.RPKM[which(rowCheck!=0),]
	if(length(unique(non0.log2.RPKM[,1]) ) != nrow(non0.log2.RPKM) ){
		print(cancer)
		break
	}
	num.log2.RPKM = non0.log2.RPKM[, -1]
	rownames(num.log2.RPKM) = non0.log2.RPKM[, 1]
	
	expr.mat = as.matrix(num.log2.RPKM)
	gsva(expr.mat, pathway.mdb.genes.list, method="ssgsea") -> es.mat
	save(es.mat, file=paste(cancer,".ssgsea.result.RData", sep=""))	
}

