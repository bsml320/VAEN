setwd("/path/to/VAEN/result.EN/dr.CCLE/")

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")
to
#####################################################################################
load("../../DATA/TCGA.ss.mat.RData")
#####################################################################################
anno = read.csv("../../DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))
#####################################################################################

solid.sample.size = solid.in_sample_R2.mat = solid.avg_CV_R2.mat = solid.F1_R2.mat = matrix(0, nrow=100, ncol=length(drugs))
colnames(solid.sample.size) = colnames(solid.in_sample_R2.mat) = colnames(solid.avg_CV_R2.mat) = colnames(solid.F1_R2.mat) = drugs

solid.mat = c()
for(ksigmoid in 1:100){
	load(paste("01S/", ksigmoid, ".CCLE.model.list.S.RData", sep=""))
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		model.list[[ drugs[kdrug] ]] -> res.list
		if(length(res.list) == 0)next
		fit <- res.list$model
		Ys = res.list$Ys
		which(Ys[,1]!=-9) -> ii
		Ys = Ys[ii, ]
		if(sd(Ys[,2]) == 0)next
		solid.mat = rbind(solid.mat, c(ksigmoid, res.list$model_summary, cor(Ys[,1], Ys[,2])) )
		
		#### way 4, PCC
		recall    = cor(Ys[,1], Ys[,2])
		precision = as.numeric(res.list$model_summary[5])
		
		solid.sample.size[ksigmoid, kdrug]      = nrow(Ys)
		solid.in_sample_R2.mat[ksigmoid, kdrug] = recall
		solid.avg_CV_R2.mat[ksigmoid, kdrug]    = precision
		solid.F1_R2.mat[ksigmoid, kdrug]        = as.numeric(res.list$model_summary[7])
	}
	cat(ksigmoid, ".", sep="")
}

save(solid.mat, solid.sample.size, solid.in_sample_R2.mat, solid.avg_CV_R2.mat, solid.F1_R2.mat, file="CCLE.S.info.RData")

############################################################################################

pdf("CCLE.S.ROC.pdf", width=5, height=5)
for(k in 1:length(drugs)){
	drug = drugs[k]
	plot(x=solid.in_sample_R2.mat[,k], y=solid.avg_CV_R2.mat[,k], main=drugs[k], xlab="Self in_sample PCC", ylab="avg PCC (in_sample)", col=rep("blue",200), pch=20, cex=.6 )
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	idx = tmp[1:10,1]
	points(solid.in_sample_R2.mat[idx,k], solid.avg_CV_R2.mat[idx,k], pch=4, col="red")
}
dev.off()

############################################################################################
############################################################################################

TCGA.pred.mat = c()
solid.model_summary = c()
holdout.R2 = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),] ### avg_CV_R2
	holdout.R2 = rbind(holdout.R2, c(drug, tmp[1,4]) )
	
	best.index = tmp[1,1]
	
	load( paste("01S/", best.index ,".CCLE.model.list.S.RData", sep="") )
	model.list[[ drug ]] -> res.list
	fit <- res.list$model
	
	TCGA.pred = read.table(paste("../../result/", best.index, ".TCGA.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
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
write.table(TCGA.pred.mat, file="VAEN_CCLE.S.pred_TCGA.txt", quote=F, sep="\t", row.names=FALSE)

############################################################################################

CCLE.PCC = c()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	pred.mat = c()
	best.index = tmp[1,1]
	
	load( paste("01S/", best.index ,".CCLE.model.list.S.RData", sep="") )
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
write.table(self.prediction.mat, file="VAEN_CCLE.S.pred_CCLE.txt", quote=F, sep="\t", row.names=FALSE)
############################################################################################

CCLE.pred.full.mat = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("01S/", best.index ,".CCLE.model.list.S.RData", sep="") )
	model.list[[ drug ]] -> res.list
	
	CCLE.latent = read.table(paste("../../result/", best.index, ".CCLE.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
	CCLE.latent.data = CCLE.latent[,-1]
	fit <- res.list$model
	CCLE.probabilities = predict(fit, as.matrix(CCLE.latent.data), s = 'lambda.min')
	
	CCLE.pred.full.mat = cbind(CCLE.pred.full.mat, CCLE.probabilities)
	cat(drug, ".", sep="")
}

self.pred.mat = cbind(CELLINE=CCLE.latent[,1], CCLE.pred.full.mat)
colnames(self.pred.mat) = c("CELLINE", drugs)
write.table(self.pred.mat, file=paste("VAEN_CCLE.S.pred_CCLE.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################################################
