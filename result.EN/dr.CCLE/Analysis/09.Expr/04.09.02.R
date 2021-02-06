setwd("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/04-mix/09.Expr/")

###############################################################################

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

###############################################################################
cur.cancer = "BRCA"

original.TCGA.RPKM = read.delim(paste("/data/mshao/TCGA/",cur.cancer,"/HiSeqV2", sep=""), as.is=T)
apply(original.TCGA.RPKM[,-1],1,sum) -> rowCheck
non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
apply(original.TCGA.RPKM[,-1],1,var) -> rowCheck
non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
TCGA.ss = gsub("\\.", "-", colnames(non0.TCGA.RPKM))
colnames(non0.TCGA.RPKM) = TCGA.ss; rm(TCGA.ss)
match(unique(non0.TCGA.RPKM[,1]), non0.TCGA.RPKM[,1]) -> ii ### keep unique genes
non0.TCGA.RPKM = non0.TCGA.RPKM[ii, ]
	
subtype = read.delim("/data1_2/jiap/projects/18-CCLE-VAE/DATA/BRCA.subtype.TCGA.txt", as.is=T, header=T, sep="\t", skip=1)
subtype[which(subtype$PAM50.mRNA=="Basal-like"),1] -> g1
subtype[which(subtype$PAM50.mRNA=="HER2-enriched"),1] -> g2
subtype[which(subtype$PAM50.mRNA=="Luminal A"),1] -> g3
subtype[which(subtype$PAM50.mRNA=="Luminal B"),1] -> g4
subtype.list = list()
subtype.list[["BRCA-Basal"]] = g1
subtype.list[["BRCA-HER2"]] = g2
subtype.list[["BRCA-LumA"]] = g3
subtype.list[["BRCA-LumB"]] = g4

for(k in 1:length(subtype.list)){
	gene.drug.expr.mat = c()
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cur.cancer), ]
	fixed.ss = intersect( substr(colnames(non0.TCGA.RPKM), 1,12), substr(blca.ccle[,1],1,12))
	fixed.ss = intersect(fixed.ss, subtype.ss)
	
	blca.ccle = blca.ccle[match(fixed.ss, substr(blca.ccle[,1], 1,12) ),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, scale) -> new.ccle
	TCGA.RPKM = non0.TCGA.RPKM[, match(fixed.ss, substr(colnames(non0.TCGA.RPKM), 1,12 )   )]
	rownames(TCGA.RPKM) = non0.TCGA.RPKM[,1]
	apply(TCGA.RPKM, 1, sd) -> rowCheck
	TCGA.RPKM = TCGA.RPKM[rowCheck!=0, ]
	
	apply(TCGA.RPKM,1,mean) -> rowMean
	TCGA.RPKM = TCGA.RPKM[rowMean > 1, ]
	
	### unique in W5
	t(apply(TCGA.RPKM,1,scale)) -> scaled.TCGA.RPKM
	dimnames(scaled.TCGA.RPKM) = dimnames(TCGA.RPKM)
	
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
		
		gene.drug.expr.mat = rbind(gene.drug.expr.mat, cbind( names(subtype.list)[k], rownames(TCGA.RPKM)[kgene], drugs, ps, tvalue, betas, stds, mean(gene.expr[which(X2==0)]), mean(gene.expr[which(X2==1)]), mean(gene.expr[which(X2==2)]) ))
		
		if(kgene %% 100 == 0)cat(kgene, ".", sep="")
	}
	cat(cur.cancer, ".", sep="")
	write.table(gene.drug.expr.mat, file=paste("W5.",names(subtype.list)[k],".Expr.txt", sep=""), quote=F, row.names=F, sep="\t")
}

###############################################################################

cur.cancer = "LGG"	
original.TCGA.RPKM = read.delim(paste("/data/mshao/TCGA/",cur.cancer,"/HiSeqV2", sep=""), as.is=T)
apply(original.TCGA.RPKM[,-1],1,sum) -> rowCheck
non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
apply(original.TCGA.RPKM[,-1],1,var) -> rowCheck
non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
TCGA.ss = gsub("\\.", "-", colnames(non0.TCGA.RPKM))
colnames(non0.TCGA.RPKM) = TCGA.ss; rm(TCGA.ss)
match(unique(non0.TCGA.RPKM[,1]), non0.TCGA.RPKM[,1]) -> ii ### keep unique genes
non0.TCGA.RPKM = non0.TCGA.RPKM[ii, ]
	
subtype = read.delim("/data1_2/jiap/projects/18-CCLE-VAE/DATA/LGG.subtype.TCGA.txt", as.is=T, header=T)
subtype[which(subtype[,9]=="coc1"),1] -> g1
subtype[which(subtype[,9]=="coc2"),1] -> g2
subtype[which(subtype[,9]=="coc3"),1] -> g3
subtype.list = list()
subtype.list[["LGG-coc1"]] = g1
subtype.list[["LGG-coc2"]] = g2
subtype.list[["LGG-coc3"]] = g3

for(k in 1:length(subtype.list)){
	gene.drug.expr.mat = c()
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cur.cancer), ]
	fixed.ss = intersect( substr(colnames(non0.TCGA.RPKM), 1,12), substr(blca.ccle[,1],1,12))
	fixed.ss = intersect(fixed.ss, subtype.ss)
	
	blca.ccle = blca.ccle[match(fixed.ss, substr(blca.ccle[,1], 1,12) ),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, scale) -> new.ccle
	TCGA.RPKM = non0.TCGA.RPKM[, match(fixed.ss, substr(colnames(non0.TCGA.RPKM), 1,12 )   )]
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
		
		gene.drug.expr.mat = rbind(gene.drug.expr.mat, cbind( names(subtype.list)[k], rownames(TCGA.RPKM)[kgene], drugs, ps, tvalue, betas, stds, mean(gene.expr[which(X2==0)]), mean(gene.expr[which(X2==1)]), mean(gene.expr[which(X2==2)]) ))
		
		if(kgene %% 100 == 0)cat(kgene, ".", sep="")
	}
	cat(cur.cancer, ".", sep="")
	write.table(gene.drug.expr.mat, file=paste("W5.",names(subtype.list)[k],".Expr.txt", sep=""), quote=F, row.names=F, sep="\t")
}

###############################################################################
cur.cancer = "THCA"

original.TCGA.RPKM = read.delim(paste("/data/mshao/TCGA/",cur.cancer,"/HiSeqV2", sep=""), as.is=T)
apply(original.TCGA.RPKM[,-1],1,sum) -> rowCheck
non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
apply(original.TCGA.RPKM[,-1],1,var) -> rowCheck
non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
TCGA.ss = gsub("\\.", "-", colnames(non0.TCGA.RPKM))
colnames(non0.TCGA.RPKM) = TCGA.ss; rm(TCGA.ss)
match(unique(non0.TCGA.RPKM[,1]), non0.TCGA.RPKM[,1]) -> ii ### keep unique genes
non0.TCGA.RPKM = non0.TCGA.RPKM[ii, ]
	
subtype = read.delim("/data1_2/jiap/projects/18-CCLE-VAE/DATA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
ss = substr(subtype[,1], 1, 15)
ss[which(subtype[,2]=="Braf-like")] -> g1
ss[which(subtype[,2]=="Ras-like")] -> g2
subtype.list = list()
subtype.list[["THCA-Braf-like"]] = g1
subtype.list[["THCA-Ras-like"]] = g2

for(k in 1:length(subtype.list)){
	gene.drug.expr.mat = c()
	cancer = names(subtype.list)[k]
	subtype.ss = substr(subtype.list[[k]], 1, 12)
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cur.cancer), ]
	fixed.ss = intersect( substr(colnames(non0.TCGA.RPKM), 1,12), substr(blca.ccle[,1],1,12))
	fixed.ss = intersect(fixed.ss, subtype.ss)
	
	blca.ccle = blca.ccle[match(fixed.ss, substr(blca.ccle[,1], 1,12) ),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, scale) -> new.ccle
	TCGA.RPKM = non0.TCGA.RPKM[, match(fixed.ss, substr(colnames(non0.TCGA.RPKM), 1,12 )   )]
	rownames(TCGA.RPKM) = non0.TCGA.RPKM[,1]
	apply(TCGA.RPKM, 1, sd) -> rowCheck
	TCGA.RPKM = TCGA.RPKM[rowCheck!=0, ]
	
	apply(TCGA.RPKM,1,mean) -> rowMean
	TCGA.RPKM = TCGA.RPKM[rowMean > 1, ]
	
	### unique in W5
	t(apply(TCGA.RPKM,1,scale)) -> scaled.TCGA.RPKM
	dimnames(scaled.TCGA.RPKM) = dimnames(TCGA.RPKM)
	
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
		
		gene.drug.expr.mat = rbind(gene.drug.expr.mat, cbind( names(subtype.list)[k], rownames(TCGA.RPKM)[kgene], drugs, ps, tvalue, betas, stds, mean(gene.expr[which(X2==0)]), mean(gene.expr[which(X2==1)]), mean(gene.expr[which(X2==2)]) ))
		
		if(kgene %% 100 == 0)cat(kgene, ".", sep="")
	}
	cat(cur.cancer, ".", sep="")
	write.table(gene.drug.expr.mat, file=paste("W5.",names(subtype.list)[k],".Expr.txt", sep=""), quote=F, row.names=F, sep="\t")
}
