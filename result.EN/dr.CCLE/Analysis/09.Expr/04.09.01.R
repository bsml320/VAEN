### changes in W5: not only Y needs to be scaled, X needs also to be scaled across samples, such that t-values across cancer types are comparable
### About X: no difference whether scale or not -- the rank remains unchanged
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("Must provide start and end\n", call.=FALSE)
} else if (length(args)==1) {
	cur.cancer = args[1]
} 
print( cur.cancer )

setwd("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/04-mix/09.Expr/")

###########################################################################################################

drug.ccle = read.table(file="/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/01/MIX-F1-W5-PCC/MIX-F1-W5-PCC.best.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, blca.ccle)
}
drug.ccle = cancer.drug.ccle

###########################################################################################################

	gene.drug.expr.mat = c()
	
	original.TCGA.RPKM = read.delim(paste("/data/mshao/TCGA/",cur.cancer,"/HiSeqV2", sep=""), as.is=T)
	
	apply(original.TCGA.RPKM[,-1],1,sum) -> rowCheck
	non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
	
	apply(original.TCGA.RPKM[,-1],1,var) -> rowCheck
	non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
	
	TCGA.ss = gsub("\\.", "-", colnames(non0.TCGA.RPKM))
	colnames(non0.TCGA.RPKM) = TCGA.ss; rm(TCGA.ss)
	
	match(unique(non0.TCGA.RPKM[,1]), non0.TCGA.RPKM[,1]) -> ii ### keep unique genes
	non0.TCGA.RPKM = non0.TCGA.RPKM[ii, ]
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cur.cancer), ]
	
	fixed.ss = intersect(colnames(non0.TCGA.RPKM), blca.ccle[,1])
	
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, scale) -> new.ccle  ### unique in W5
	TCGA.RPKM = non0.TCGA.RPKM[, match(fixed.ss, colnames(non0.TCGA.RPKM))]
	rownames(TCGA.RPKM) = non0.TCGA.RPKM[,1]
	apply(TCGA.RPKM, 1, sd) -> rowCheck
	TCGA.RPKM = TCGA.RPKM[rowCheck!=0, ]
	
	apply(TCGA.RPKM,1,mean) -> rowMean
	TCGA.RPKM = TCGA.RPKM[rowMean > 1, ]
	
	for(kgene in 1:nrow(TCGA.RPKM)){
		t(TCGA.RPKM[kgene, ]) -> gene.expr
		
		X2 = rep(1, length(gene.expr))
		X2[which(gene.expr < quantile(gene.expr, probs=.25))] = 0
		X2[which(gene.expr > quantile(gene.expr, probs=.75))] = 2
		if(length(unique(X2)) < 3)next
		
		apply(new.ccle, 2, function(Y){
				summary(glm(Y ~ X2)) -> sfit
				coef(sfit)[2, 1:4]
		}) -> ps.twosided
		
		betas  = ps.twosided[1,]
		stds   = ps.twosided[2,]
		tvalue = ps.twosided[3,]
		ps     = ps.twosided[4,]
		
		gene.drug.expr.mat = rbind(gene.drug.expr.mat, cbind(cur.cancer, rownames(TCGA.RPKM)[kgene], drugs, ps, tvalue, betas, stds, mean(gene.expr[which(X2==0)]), mean(gene.expr[which(X2==1)]), mean(gene.expr[which(X2==2)]) ))
		
		if(kgene %% 100 == 0)cat(kgene, ".", sep="")
	}
	cat(cur.cancer, ".", sep="")
	write.table(gene.drug.expr.mat, file=paste("W5.",cur.cancer,".Expr.txt", sep=""), quote=F, row.names=F, sep="\t")


