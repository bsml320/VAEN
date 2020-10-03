setwd("/work/RANK.Sigmoid/result.EN/dr.CCLE/01/F1-W5-PCC/")

library(MASS)
library("magrittr")
library("glmnet")
library("modEvA")
library(vegan)

calc_R2 <- function(y, y_pred) {
	tss <- sum(y**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

calc_MSE <- function(y, y_pred) {
	mean( (y - y_pred)^2 )
}

#####################################################################################
load("/work/data/TCGA.ss.mat.RData")
#####################################################################################
anno = read.csv("/work/data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))
#####################################################################################

all.sample.size = all.in_sample_R2.mat = all.avg_CV_R2.mat = all.F1_R2.mat = matrix(0, nrow=100, ncol=length(drugs))
colnames(all.sample.size) = colnames(all.in_sample_R2.mat) = colnames(all.avg_CV_R2.mat) = colnames(all.F1_R2.mat) = drugs

all.mat = c()
for(ksigmoid in 1:100){
	fn = paste("/work/RANK.Sigmoid/result.EN/dr.CCLE/01/", ksigmoid, ".CCLE.model.list.RData", sep="")
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
		
		all.in_sample_R2.mat[ksigmoid, kdrug] = cor(Ys[,1], Ys[,2])
		all.avg_CV_R2.mat[ksigmoid, kdrug]    = as.numeric(res.list$model_summary[5])
		all.F1_R2.mat[ksigmoid, kdrug]        = as.numeric(res.list$model_summary[7])
		
	}
	cat(ksigmoid, ".", sep="")
}

############################################################################################
############################################################################################
### best

TCGA.pred.mat = c()
all.model_summary = c()
holdout.R2 = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),] ### avg_CV_R2
	holdout.R2 = rbind(holdout.R2, c(drug, tmp[1,4]) )
	
	best.index = tmp[1,1]
	
	load( paste("/work/RANK.Sigmoid/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
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

CCLE.PCC = c()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	pred.mat = c()
	best.index = tmp[1,1]
	
	load( paste("/work/RANK.Sigmoid/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	Ys = res.list$Ys
		
	if(kdrug == 1){
		self.prediction.mat = Ys
	} else {
		self.prediction.mat = cbind(self.prediction.mat, Ys[,2])
	}
}

self.prediction.mat[,1] = rownames(Ys)
colnames(self.prediction.mat) = c("CELLLINE", drugs)
write.table(self.prediction.mat, file="F1-W5-PCC.best.pred_CCLE.txt", quote=F, sep="\t", row.names=FALSE)
############################################################################################

CCLE.pred.full.mat = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	pred.mat = c()
	best.index = tmp[1,1]
	
	load(                  paste("/work/RANK.Sigmoid/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	
	CCLE.latent = read.table(paste("/work/RANK.Sigmoid/result/", best.index, ".CCLE.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
	CCLE.latent.data = CCLE.latent[,-1]
	fit <- res.list$model
	CCLE.probabilities = predict(fit, as.matrix(CCLE.latent.data), s = 'lambda.min')
	pred.mat = cbind(pred.mat, scale(CCLE.probabilities))
	
	CCLE.pred.full.mat = cbind(CCLE.pred.full.mat, pred.mat)
	cat(drug, ".", sep="")
}

self.pred.mat = cbind(CELLINE=CCLE.latent[,1], CCLE.pred.full.mat)
colnames(self.pred.mat) = c("CELLINE", drugs)
write.table(self.pred.mat, file=paste("F1-W5-PCC.best.pred_CCLE.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################################################
