setwd("/work/RANK.Sigmoid/result.EN/dr.GDSC/01/F1-W5-PCC/")

library(MASS)
library("magrittr")
library("glmnet")
library("modEvA")
library(vegan)

#####################################################################################
load("/work/data/TCGA.ss.mat.RData")
#####################################################################################
anno = read.delim("/work/data/GDSC/v17.3_fitted_dose_response.txt", as.is=T)
drugs = sort(unique(anno$DRUG_NAME))
#####################################################################################

all.sample.size = all.in_sample_R2.mat = all.avg_CV_R2.mat = all.F1_R2.mat = matrix(0, nrow=100, ncol=length(drugs))
colnames(all.sample.size) = colnames(all.in_sample_R2.mat) = colnames(all.avg_CV_R2.mat) = colnames(all.F1_R2.mat) = drugs

all.mat = c()
for(ksigmoid in 1:100){
	fn = paste("/work/RANK.Sigmoid/result.EN/dr.GDSC/01/", ksigmoid, ".GDSC.model.list.RData", sep="")
	if(!file.exists(fn))next
	load(fn)
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		model.list[[ drugs[kdrug] ]] -> res.list
		if(length(res.list) == 0)next
		fit <- res.list$model
		Ys = res.list$Ys
		which(Ys[,1]!=-9) -> ii
		Ys = Ys[ii, ]
		if(sd(Ys[,2]) == 0)next
		all.mat = rbind(all.mat, c(ksigmoid, res.list$model_summary, cor(Ys[,1], Ys[,2])) )
		
		#### way 4, PCC
		recall    = cor(Ys[,1], Ys[,2])
		precision = as.numeric(res.list$model_summary[5])
		
		all.sample.size[ksigmoid, kdrug]      = nrow(Ys)
		all.in_sample_R2.mat[ksigmoid, kdrug] = recall
		all.avg_CV_R2.mat[ksigmoid, kdrug]    = precision
		#all.F1_R2.mat[ksigmoid, kdrug]        = 2 * precision * recall/(precision + recall)
		all.F1_R2.mat[ksigmoid, kdrug]        = as.numeric(res.list$model_summary[7])
	}
	cat(ksigmoid, ".", sep="")
}

#####################################################################################
### best

TCGA.pred.mat = c()
all.model_summary = c()
holdout.R2 = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	best.index = tmp[1,1]
	
	load( paste("/work/RANK.Sigmoid/result.EN/dr.GDSC/01/", best.index ,".GDSC.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	fit <- res.list$model
	
	TCGA.pred = read.table(paste("/work/RANK.Sigmoid/result/", best.index, ".TCGA.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
	TCGA.test.data = TCGA.pred[,-1]
	TCGA.probabilities = predict(fit, as.matrix(TCGA.test.data), s = 'lambda.min')
	
	TCGA.pred.mat = cbind(TCGA.pred.mat, TCGA.probabilities)
	cat("...", drug, ".", sep="")
}

TCGA.pred.mat = cbind(TCGA.pred[,1], "A", TCGA.pred.mat)
gsub("\\.", "-", TCGA.pred.mat[,1]) -> ss
TCGA.pred.mat[,1] = ss
match(TCGA.pred.mat[,1], TCGA.ss.mat[,1]) -> ii
TCGA.pred.mat[,2] = TCGA.ss.mat[ii, 2]
colnames(TCGA.pred.mat) = c("Sample", "Cancer", drugs)
write.table(TCGA.pred.mat, file="F1-W5-PCC.best.pred_TCGA.txt", quote=F, sep="\t", row.names=FALSE)
############################################################################################

GDSC.pred.full.mat = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	best.index = tmp[1,1]
	load(                    paste("/work/RANK.Sigmoid/result.EN/dr.GDSC/01/", best.index ,".GDSC.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
		
	GDSC.latent = read.table(paste("/work/RANK.Sigmoid/result/", best.index, ".CCLE.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
	GDSC.latent.data = GDSC.latent[,-1]
	fit <- res.list$model
	GDSC.probabilities = predict(fit, as.matrix(GDSC.latent.data), s = 'lambda.min')
	
	GDSC.pred.full.mat = cbind(GDSC.pred.full.mat, GDSC.probabilities)
	cat(drug, ".", sep="")
}

self.pred.mat = cbind(CELLINE=GDSC.latent[,1], GDSC.pred.full.mat)
colnames(self.pred.mat) = c("CELLINE", drugs)
write.table(self.pred.mat, file=paste("F1-W5-PCC.best.pred_GDSC.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################################################
############################################################################################

calc_R2 <- function(y, y_pred) {
	tss <- sum(y**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

GDSC.R2 = c()
GDSC.PCC = c()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("/work/RANK.Sigmoid/result.EN/dr.GDSC/01/", best.index ,".GDSC.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	Ys = res.list$Ys
	which(Ys[,1]!=-9) -> ii
	calc_R2(Ys[ii, 1], Ys[ii, 2]) -> r2
	GDSC.R2 = rbind(GDSC.R2, c(drug, r2) )
	cor(Ys[ii, 1], Ys[ii, 2]) -> pcc
	GDSC.PCC = rbind(GDSC.PCC, c(drug, pcc) )
	
	if(kdrug == 1){
		self.prediction.mat = Ys
	} else {
		self.prediction.mat = cbind(self.prediction.mat, Ys[,2])
	}
	cat(kdrug,".",sep="")
}
write.table(GDSC.R2,  file="F1-W5-PCC.best.GDSC.R2.txt", quote=F, sep="\t", row.names=FALSE)
write.table(GDSC.PCC, file="F1-W5-PCC.best.GDSC.PCC.txt", quote=F, sep="\t", row.names=FALSE)

self.prediction.mat[,1] = rownames(Ys)
colnames(self.prediction.mat) = c("CELLLINE", drugs)
write.table(self.prediction.mat, file="F1-W5-PCC.best.pred_GDSC.txt", quote=F, sep="\t", row.names=FALSE)

############################################################################################
