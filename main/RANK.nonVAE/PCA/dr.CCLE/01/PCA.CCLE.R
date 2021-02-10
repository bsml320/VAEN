setwd("/path/to/VAEN/main/RANK.nonVAE/PCA/dr.CCLE/01")

dat = read.table("/path/to/VAEN/main/V15.CCLE.4VAE.RANK.tsv", header=T, as.is=T)
res.pca <- prcomp(dat, scale = TRUE)

PPs = res.pca$x[, 1:100]

TCGA = read.table("/path/to/VAEN/main/V15.TCGA.4VAE.RANK.tsv", header=T, as.is=T)
ind.sup.coord <- predict(res.pca, newdata = TCGA)

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")

source("/path/to/VAEN/code/nested_EN.R")

library(parallel)
parallel.main = function(kk, train.data, Y, n_folds=10, n_train_test_folds=5, seed=NA, alpha=0.5, null_testing=FALSE, drug){
	main(train.data, Y, n_folds=n_folds, n_train_test_folds=n_train_test_folds, seed=seed, alpha=alpha, null_testing=null_testing, drug=drug) -> res.list
	res.list
}

#####################################################################################
load("/path/to/VAEN/DATA/TCGA.ss.mat.RData")
#####################################################################################

	########### Prediction

	original.ss.PP = rownames(PPs)
	sapply(original.ss.PP, function(x){
		new.u = u = strsplit(x, split="\\.")[[1]][1]
		if(grepl("^X", u)){
			substr(u, 2, nchar(u)) -> new.u
		}
		new.u
	}) -> ss.PP
	names(ss.PP) = NULL

	########### original drug data

	anno = read.csv("/path/to/VAEN/DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
	drugs = sort(unique(anno$Compound))

	model.list = list()
	for(k in 1:length(drugs)){
		drug = drugs[k]
		cat(drugs[k], " ======== start\n", sep="")
		
		anno.1 = anno[which(anno$Compound==drugs[k] & anno[, "ActArea"] != -9 ),]
		intersect(anno.1[,1], ss.PP) -> shared.samples
		match(shared.samples, anno.1[,1]) -> ii.Y
		anno.2 = anno.1[ii.Y, ]
		Y = anno.2[, "ActArea"]
		
		match(shared.samples, ss.PP) -> ii
		train.data = PPs[ii,]
		
		tmp.list = list()
		mclapply(1:10, parallel.main, train.data, Y, n_folds=10, n_train_test_folds=5, seed=NA, alpha=0.5, null_testing=FALSE, drug=drugs[k], mc.cores=10) -> test
		for(kk in 1:10){
			res.list = test[[kk]]
			model = res.list$model
			beta = coef(model, model$lambda.min)
			n = sum(beta!=0)
			cat( "PCC = ", cor(Y, res.list$self_pred[,1]), "; n = ", n, "\n", sep="")
			tmp.list[[kk]] = res.list
		}
		
		unlist(lapply(tmp.list, function(u)as.numeric(u$model_summary[5]))) -> cv_R2_avg
		which.max(cv_R2_avg) -> idx
		
		res.list = tmp.list[[idx]]
		cat("selected ", idx, ", ", sep="")
		
		Ys = matrix(-9, nrow=nrow(PPs), ncol=2)
		rownames(Ys) = rownames(PPs)
		Ys[ match(shared.samples, ss.PP) ,1] = Y
		Ys[ match(shared.samples, ss.PP) ,2] = res.list$self_pred[,1]
		res.list$Ys = Ys
		model.list[[ drugs[k] ]] = res.list
		
		cat(drugs[k]," end \n", sep="")
	}

	########### match
	save(model.list, file=paste("PCA.CCLE.model.list.RData", sep=""))

#####################################################################################
anno = read.csv("/path/to/VAEN/DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))
#####################################################################################

	info = c()
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		model.list[[ drugs[kdrug] ]] -> res.list
		if(length(res.list) == 0)next
		fit <- res.list$model
		Ys = res.list$Ys
		which(Ys[,1]!=-9) -> ii
		Ys = Ys[ii, ]
		if(sd(Ys[,2]) == 0)next
		
		#### way 4, PCC
		recall    = cor(Ys[,1], Ys[,2])
		precision = as.numeric(res.list$model_summary[5])
		F1_R2 = 2 * precision * recall/(precision + recall)
		info = rbind(info, c(drug, F1_R2, recall, precision) )
	}

colnames(info) = c("Drug", "F1", "Self_PCC", "avg_CV_R2")

write.table(info, file="NOPEER.RANK.PCA.best.txt", row.names=F, col.names=T, quote=F, sep="\t")

#####################################################################################
#####################################################################################

TCGA.PPs = ind.sup.coord[, 1:100]
CCLE.latent.data = res.pca$x[, 1:100]
TCGA.pred.mat = c()
CCLE.pred.full.mat = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	
	model.list[[ drug ]] -> res.list
	fit <- res.list$model
	
	
	TCGA.probabilities = predict(fit, as.matrix(TCGA.PPs), s = 'lambda.min')
	TCGA.pred.mat = cbind(TCGA.pred.mat, TCGA.probabilities)
	
	
	CCLE.probabilities = predict(fit, as.matrix(CCLE.latent.data), s = 'lambda.min')
	CCLE.pred.full.mat = cbind(CCLE.pred.full.mat, CCLE.probabilities)
	
	Ys = res.list$Ys
	if(k == 1){
		self.prediction.mat = Ys
	} else {
		self.prediction.mat = cbind(self.prediction.mat, Ys[,2])
	}
	
	cat("...", drug, ".", sep="")
}

TCGA.pred.mat = cbind(rownames(TCGA.PPs), "A", TCGA.pred.mat)
gsub("\\.", "-", TCGA.pred.mat[,1]) -> ss
TCGA.pred.mat[,1] = ss
match(TCGA.pred.mat[,1], TCGA.ss.mat[,1]) -> ii
TCGA.pred.mat[,2] = TCGA.ss.mat[ii, 2]
colnames(TCGA.pred.mat) = c("Sample", "Cancer", drugs)
write.table(TCGA.pred.mat, file="PCA.pred_TCGA.txt", quote=F, sep="\t", row.names=FALSE)

self.prediction.mat[,1] = rownames(Ys)
colnames(self.prediction.mat) = c("CELLLINE", drugs)
write.table(self.prediction.mat, file="PCA.pred_CCLE.txt", quote=F, sep="\t", row.names=FALSE)

self.pred.mat = cbind(CELLINE=rownames(CCLE.latent.data), CCLE.pred.full.mat)
colnames(self.pred.mat) = c("CELLINE", drugs)
write.table(self.pred.mat, file=paste("PCA.pred_CCLE.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

