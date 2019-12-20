### Prepare data for Figure 3BCD
### This script will generate two files: CCLE.ActArea.lineage.RData and PLX4720.RData

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")

setwd("/work/Figures/Figure3/")
load("/work/data/TCGA.ss.mat.RData")

####################################### obsd
load(paste("/work/result.EN/dr.CCLE/01/1.CCLE.model.list.RData", sep=""))
ccle.model.list = model.list
ccle.obsd = c()
for(k in 1:length(ccle.model.list)){
	res.list = ccle.model.list[[k]]
	drug = names(ccle.model.list)[k]
	Ys = res.list$Ys
	
	ccle.obsd = cbind(ccle.obsd, Ys[,1])
}
colnames(ccle.obsd) = names(ccle.model.list)
#######################################

ccle.pred      = read.table(file="/work/result.EN/dr.CCLE/01/F1-W5-PCC/F1-W5-PCC.avgtop10.pred_CCLE.txt", header=T, as.is=T, sep="\t")
ccle.pred.full = read.table(file="/work/result.EN/dr.CCLE/01/F1-W5-PCC/F1-W5-PCC.avgtop10.pred_CCLE.full.txt", header=T, as.is=T, sep="\t")

colnames(ccle.pred) -> ss
ss[which(ss=="X17.AAG")] = "17-AAG"
gsub("\\.", "-", ss) -> ss; colnames(ccle.pred) = ss

colnames(ccle.pred.full) -> ss
ss[which(ss=="X17.AAG")] = "17-AAG"
gsub("\\.", "-", ss) -> ss; colnames(ccle.pred.full) = ss

#####################################################################################
########### Prediction
original.ss.PP = rownames(ccle.obsd)
sapply(original.ss.PP, function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
names(tt) = NULL

#####################################################
anno = read.csv("/work/data/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))

full.drug.by.tissue.mat = pred.drug.by.tissue.mat = obsd.drug.by.tissue.mat = matrix(1, nrow=length(drugs), ncol=length(unique(tt)))
rownames(obsd.drug.by.tissue.mat) = rownames(pred.drug.by.tissue.mat) = rownames(full.drug.by.tissue.mat) = drugs
colnames(obsd.drug.by.tissue.mat) = colnames(pred.drug.by.tissue.mat) = colnames(full.drug.by.tissue.mat) = sort(unique(tt))

for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	Ys = cbind(ccle.obsd[, drug], ccle.pred[, drug], ccle.pred.full[, drug])
	
	sapply(rownames(Ys), function(x){
		strsplit(x, split="\\.")[[1]][1] -> u
		strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
	}) -> tt
	names(tt) = NULL
	
	for(kt in 1:ncol(obsd.drug.by.tissue.mat)){
		x.tissue = colnames(obsd.drug.by.tissue.mat)[kt]
		
		X = ifelse(tt %in% x.tissue, 1, 0)
		
		which(Ys[,1]!=-9) -> ii
		Y1 = Ys[ii, 1]
		X1 = X[ii]
		sum(X1==1) -> check
		label = ifelse(mean(Y1[which(X1==1)]) > mean(Y1[which(X1==0)]), 1, -1)
		if(check > 5)obsd.drug.by.tissue.mat[kdrug, kt] = label * t.test(Y1[which(X1==1)], Y1[which(X1==0)])$p.value
		
		which(Ys[,1]!=-9) -> ii
		Y2 = Ys[ii, 2]
		X2 = X[ii]
		sum(X2==1) -> check
		label = ifelse(mean(Y2[which(X2==1)]) > mean(Y2[which(X2==0)]), 1, -1)
		if(check > 5)pred.drug.by.tissue.mat[kdrug, kt] = label * t.test(Y2[which(X2==1)], Y2[which(X2==0)])$p.value
		
		Y3 = Ys[,3]
		label = ifelse(mean(Y3[which(X==1)]) > mean(Y3[which(X==0)]), 1, -1)
		full.drug.by.tissue.mat[kdrug, kt] = label * t.test(Y3[which(X==1)], Y3[which(X==0)])$p.value
	}
}


apply(obsd.drug.by.tissue.mat, 2, var) -> obsd.check
apply(pred.drug.by.tissue.mat, 2, var) -> pred.check

print(which(obsd.check==0))
print(which(pred.check==0))

obsd.drug.by.tissue.mat = obsd.drug.by.tissue.mat[, which(obsd.check != 0)]
pred.drug.by.tissue.mat = pred.drug.by.tissue.mat[, which(obsd.check != 0)]
full.drug.by.tissue.mat = full.drug.by.tissue.mat[, which(obsd.check != 0)]

########################################################################################################################
predicted.dr = read.table(paste("/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", sep=""), header=T, as.is=T)
drugs = colnames(predicted.dr)[c(-1, -2)]
cancer.types = unique(predicted.dr[,2])
sample.type = substr(predicted.dr[,1], 14, 15)

cancer.predicted.dr = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = predicted.dr[which(predicted.dr[,2] == cancer & sample.type == type.code), ]
	cancer.predicted.dr = rbind(cancer.predicted.dr, blca.ccle)
}

for(k in 3:ncol(predicted.dr)){
	if( sum(predicted.dr[,k] < 0) != 0){
		min(predicted.dr[,k]) -> shift.factor
		predicted.dr[,k] = predicted.dr[,k] + abs(shift.factor)
	}
}

##########
TCGA.drug.by.cancer.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(TCGA.drug.by.cancer.mat) = drugs
colnames(TCGA.drug.by.cancer.mat) = sort(cancer.types)

for(kdrug in 1:length(drugs)){
	Y = cancer.predicted.dr[, drugs[kdrug]]
	
	for(kt in 1:ncol(TCGA.drug.by.cancer.mat)){
		x.cancer = colnames(TCGA.drug.by.cancer.mat)[kt]
		X = ifelse(cancer.predicted.dr[,2] %in% x.cancer, 1, 0)
		if(sum(X==1) < 5)next
		label = ifelse(mean(Y[which(X==1)]) > mean(Y[which(X==0)]), 1, -1)
		TCGA.drug.by.cancer.mat[kdrug, kt] = label * t.test(Y[which(X==1)], Y[which(X==0)])$p.value
	}
}

save(obsd.drug.by.tissue.mat, pred.drug.by.tissue.mat, full.drug.by.tissue.mat, TCGA.drug.by.cancer.mat, file="CCLE.ActArea.lineage.RData")

########################################################################################################################
TCGA.sensitive.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(TCGA.sensitive.mat) = drugs
colnames(TCGA.sensitive.mat) = sort(cancer.types)
sensitive.prop = c()

for(kdrug in 1:length(drugs)){
	Y = cancer.predicted.dr[, drugs[kdrug]]
	Y > quantile(Y, probs=.95) -> sensitive.ii
	for(kt in 1:ncol(TCGA.drug.by.cancer.mat)){
		x.cancer = colnames(TCGA.drug.by.cancer.mat)[kt]
		X = cancer.predicted.dr[,2] %in% x.cancer
		if(kdrug == 18)sensitive.prop = c(sensitive.prop, sum(X & sensitive.ii)/sum(X) )
		
		if(sum(X & sensitive.ii) < 5)next
		
		table(sensitive.ii, X) -> mm
		mm[2:1, 2:1] -> new.mm
		TCGA.sensitive.mat[kdrug, kt] = fisher.test(new.mm)$p.value
	}
}
names(sensitive.prop) = colnames(TCGA.drug.by.cancer.mat)

########################################################################################################################
TCGA.resistant.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(TCGA.resistant.mat) = drugs
colnames(TCGA.resistant.mat) = sort(cancer.types)
resistant.prop = c()

for(kdrug in 1:length(drugs)){
	Y = cancer.predicted.dr[, drugs[kdrug]]
	Y < quantile(Y, probs=.05) -> resistant.ii
	for(kt in 1:ncol(TCGA.drug.by.cancer.mat)){
		x.cancer = colnames(TCGA.drug.by.cancer.mat)[kt]
		X = cancer.predicted.dr[,2] %in% x.cancer
		if(kdrug == 18)resistant.prop = c(resistant.prop, sum(X & resistant.ii)/sum(X) )
		if(sum(X & resistant.ii) < 5)next
		
		table(resistant.ii, X) -> mm
		mm[2:1, 2:1] -> new.mm
		TCGA.resistant.mat[kdrug, kt] = fisher.test(new.mm)$p.value
	}
}
names(resistant.prop) = colnames(TCGA.drug.by.cancer.mat)

PLX4720.pred.dr = cancer.predicted.dr[, 18]
save(resistant.prop, sensitive.prop, PLX4720.pred.dr, file="PLX4720.RData")

rownames(TCGA.sensitive.mat) = rownames(pred.drug.by.tissue.mat)
rownames(TCGA.resistant.mat) = rownames(pred.drug.by.tissue.mat)

save(obsd.drug.by.tissue.mat, pred.drug.by.tissue.mat, full.drug.by.tissue.mat, TCGA.drug.by.cancer.mat, TCGA.sensitive.mat, TCGA.resistant.mat, file="CCLE.ActArea.lineage.RData")

