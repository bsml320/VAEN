### model selection and the final predicted drug response
setwd("/work/result.EN/dr.CCLE/01/final/")

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")

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
	fn = paste("/work/result.EN/dr.CCLE/01/", ksigmoid, ".CCLE.model.list.RData", sep="")
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
		
		#### In our original analyses, we calculated several parameters. In the final results, we only used avg_CV_R2 for model selection.
		recall    = cor(Ys[,1], Ys[,2])
		precision = as.numeric(res.list$model_summary[5])
		
		all.sample.size[ksigmoid, kdrug]      = nrow(Ys)
		all.in_sample_R2.mat[ksigmoid, kdrug] = recall
		all.avg_CV_R2.mat[ksigmoid, kdrug]    = precision
		all.F1_R2.mat[ksigmoid, kdrug]        = 2 * precision * recall/(precision + recall)
	}
	cat(ksigmoid, ".", sep="")
}

#####################################################################################
#####################################################################################

TCGA.pred.mat = c()
all.model_summary = c()
holdout.R2 = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),] ### avg_CV_R2
	holdout.R2 = rbind(holdout.R2, c(drug, tmp[1,4]) )
	
	candi = tmp[1:10,1]
	pred.mat = c()
	for(best.index in candi){
		load( paste("/work/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
		model.list[[ drug ]] -> res.list
		model = res.list$model
		beta = coef(model, model$lambda.min)
		
		TCGA.pred = read.table(paste("/work/result/", best.index, ".TCGA.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
		TCGA.test.data = TCGA.pred[,-1]
		fit <- res.list$model
		TCGA.probabilities = predict(fit, as.matrix(TCGA.test.data), s = 'lambda.min')
		
		pred.mat = cbind(pred.mat, scale(TCGA.probabilities))
		cat(".", sep="")
	}
	apply(pred.mat, 1, mean) -> TCGA.probabilities
	cor(pred.mat) -> mm
	mm[lower.tri(mm)] = NA
	diag(mm) = NA
	cat(mean(as.vector(mm), na.rm=T))
	
	TCGA.pred.mat = cbind(TCGA.pred.mat, TCGA.probabilities)
	cat("...", drug, ".", sep="")
}

TCGA.pred.mat = cbind(TCGA.pred[,1], "A", TCGA.pred.mat)
gsub("\\.", "-", TCGA.pred.mat[,1]) -> ss
TCGA.pred.mat[,1] = ss
match(TCGA.pred.mat[,1], TCGA.ss.mat[,1]) -> ii
TCGA.pred.mat[,2] = TCGA.ss.mat[ii, 2]
colnames(TCGA.pred.mat) = c("Sample", "Cancer", drugs)
write.table(TCGA.pred.mat, file="F1-W5-PCC.avgtop10.pred_TCGA.txt", quote=F, sep="\t", row.names=FALSE)

write.table(holdout.R2, file="CCLE.best.holdout.R2.txt", quote=F, sep="\t", row.names=FALSE)

############################################################################################

pdf("F1-W5-PCC.ROC.pdf", width=5, height=5)
for(k in 1:length(drugs)){
	drug = drugs[k]
	plot(x=all.in_sample_R2.mat[,k], y=all.avg_CV_R2.mat[,k], main=drugs[k], xlab="Self in_sample PCC", ylab="avg PCC (in_sample)", col=rep("lightgreen",200), pch=20, cex=.6 )
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	idx = tmp[1:10,1]
	points(all.in_sample_R2.mat[idx,k], all.avg_CV_R2.mat[idx,k], pch=4, col="red")
}
dev.off()

############################################################################################

CCLE.PCC = c()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	pred.mat = c()
	for(best.index in tmp[1:10,1]){
		load( paste("/work/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
		model.list[[ drug ]] -> res.list
		Ys = res.list$Ys
		which(Ys[,1]!=-9) -> ii
		
		avg = mean(Ys[ii, 2])
		sd2 = sd(Ys[ii, 2])
		scaled.Y = (Ys[,2] - avg)/sd2
		scaled.Y[-ii] = -9
		
		pred.mat = cbind(pred.mat, scaled.Y)
	}
	Ys[,2] = apply(pred.mat, 1, mean) 
	which(Ys[,1]!=-9) -> ii
	cor(Ys[ii, 1], Ys[ii, 2]) -> pcc
	CCLE.PCC = rbind(CCLE.PCC, c(drug, pcc) )
	
	if(kdrug == 1){
		self.prediction.mat = Ys
	} else {
		self.prediction.mat = cbind(self.prediction.mat, Ys[,2])
	}
}

write.table(CCLE.PCC, file="F1-W5-PCC.avgtop10.CCLE.PCC.txt", quote=F, sep="\t", row.names=FALSE)
self.prediction.mat[,1] = rownames(Ys)
colnames(self.prediction.mat) = c("CELLLINE", drugs)
write.table(self.prediction.mat, file="F1-W5-PCC.avgtop10.pred_CCLE.txt", quote=F, sep="\t", row.names=FALSE)

############################################################################################

calc_R2 <- function(y, y_pred) {
	tss <- sum(y**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

CCLE.R2 = c()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("/work/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	Ys = res.list$Ys
	
	which(Ys[,1]!=-9) -> ii
	calc_R2(Ys[ii, 1], Ys[ii, 2]) -> r2
	CCLE.R2 = rbind(CCLE.R2, c(drug, r2) )
}

write.table(CCLE.R2, file="F1-W5-PCC.best.CCLE.R2.txt", quote=F, sep="\t", row.names=FALSE)

############################################################################################
############################################################################################

CCLE.pred.full.mat = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	
	pred.mat = c()
	for(best.index in tmp[1:10,1]){
		load(                  paste("/work/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
		model.list[[ drug ]] -> res.list
		
		CCLE.latent = read.table(paste("/work/result/", best.index, ".CCLE.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
		CCLE.latent.data = CCLE.latent[,-1]
		fit <- res.list$model
		CCLE.probabilities = predict(fit, as.matrix(CCLE.latent.data), s = 'lambda.min')
		pred.mat = cbind(pred.mat, scale(CCLE.probabilities))
	}
	apply(pred.mat, 1, mean) -> CCLE.probabilities
	
	CCLE.pred.full.mat = cbind(CCLE.pred.full.mat, CCLE.probabilities)
	cat(drug, ".", sep="")
}

self.pred.mat = cbind(CELLINE=CCLE.latent[,1], CCLE.pred.full.mat)
colnames(self.pred.mat) = c("CELLINE", drugs)
write.table(self.pred.mat, file=paste("F1-W5-PCC.avgtop10.pred_CCLE.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################################################
save(all.mat, all.sample.size, all.in_sample_R2.mat, all.avg_CV_R2.mat, all.F1_R2.mat,F1.all.res.mat, file="F1-W5-PCC.info.RData")
