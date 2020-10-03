#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("Must provide start and end\n", call.=FALSE)
} else if (length(args)==2) {
	start = args[1]
	end = args[2]
	scale.factor = 1
} else if (length(args)==3) {
	start = args[1]
	end = args[2]
	scale.factor = as.numeric(args[3])
} 
print( c(start, end ) )

### for parallel running, this script can be implemented as "Rscript 01.01.CCLE.EN_GLM.R 1 10" for the VAE models 1-10, 
### "Rscript 01.01.CCLE.EN_GLM.R 11 20" for the VAE models 11-20, 
### and so on



library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")


### using RANK.Sigmoid as an example; 
### this folder can be any of RANK.Sigmoid, RANK.ReLU, ZS.Sigmoid, ZS.ReLU, Z01.Sigmoid, Z01.ReLU

setwd("/work/RANK.Sigmoid/result.EN/dr.CCLE/01")
source("/work/code/nested_EN.R")


library(parallel)
parallel.main = function(kk, train.data, Y, n_folds=10, n_train_test_folds=5, seed=NA, alpha=0.5, null_testing=FALSE, drug){
	main(train.data, Y, n_folds=n_folds, n_train_test_folds=n_train_test_folds, seed=seed, alpha=alpha, null_testing=null_testing, drug=drug) -> res.list
	res.list
}

#####################################################################################
load("/work/data/TCGA.ss.mat.RData")
#####################################################################################

for(ksigmoid in start:end){
	
	cat("ksigmoid = ", ksigmoid, "\n", sep="")

	#############
	model_summary_file <- paste(ksigmoid, ".model_summary.txt", sep="")
	model_summary_cols <- c('Drug', 'alpha', 'n_snps_in_model', 'lambda_min_mse',
							  'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
							  'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
							  'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')

	#####################################################################################

	TCGA.pred = read.table(paste("//work/RANK.Sigmoid/result/", ksigmoid, ".TCGA.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
	TCGA.test.data = TCGA.pred[,-1]

	########### Prediction

	PPs = read.table(paste("/work/RANK.Sigmoid/result/",ksigmoid,".CCLE.latent.tsv", sep=""))
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

	anno = read.csv("/work/data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
	drugs = sort(unique(anno$Compound))

	self.prediction.mat = matrix(-9, nrow=length(unique(anno[,1])), ncol=length(drugs)+2)
	self.prediction.mat[,1] = unique(anno[,"Primary.Cell.Line.Name"])
	self.prediction.mat[,2] = unique(anno[,"CCLE.Cell.Line.Name"])
	colnames(self.prediction.mat) = c("CELLLINE", "Type", drugs)

	write(model_summary_cols, file = model_summary_file, ncol = 21, sep = '\t')

	TCGA.drug.response.mat = c()
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
		
		unlist(lapply(tmp.list, function(u)as.numeric(u$model_summary[4]))) -> lambda
		unlist(lapply(tmp.list, function(u)as.numeric(u$model_summary[9]))) -> self_r2
		unlist(lapply(tmp.list, function(u)as.numeric(u$model_summary[3]))) -> n
		cbind(n, lambda, cv_R2_avg, self_r2) -> x
		x[order(x[,3]), ]

		
		if(sum(is.na(cv_R2_avg)) == length(cv_R2_avg)){
			model_summary <- c(drug, 0.5, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
			write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
			TCGA.probabilities = rep(0, nrow(TCGA.test.data))
		} else {
			res.list = tmp.list[[idx]]
			cat("selected ", idx, ", ", sep="")
			write(res.list$model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
			fit <- res.list$model
			TCGA.probabilities = predict(fit, as.matrix(TCGA.test.data), s = 'lambda.min')
		}
		TCGA.drug.response.mat = cbind(TCGA.drug.response.mat, TCGA.probabilities)
		
		Ys = matrix(-9, nrow=nrow(PPs), ncol=2)
		rownames(Ys) = rownames(PPs)
		Ys[ match(shared.samples, ss.PP) ,1] = Y
		Ys[ match(shared.samples, ss.PP) ,2] = res.list$self_pred[,1]
		res.list$Ys = Ys
		model.list[[ drugs[k] ]] = res.list
		
		match(anno.2[,1], self.prediction.mat[,2]) -> pii.1
		self.prediction.mat[pii.1, k+2] = res.list$self_pred
		
		cat(drugs[k]," end \n", sep="")
	}

	########### match
	gsub("\\.", "-", TCGA.pred[,1]) -> ss; TCGA.pred[,1] = ss
	match(TCGA.pred[,1], TCGA.ss.mat[,1]) -> ii
	TCGA.drug.response.mat = cbind(TCGA.pred[,1], TCGA.ss.mat[ii, 2], TCGA.drug.response.mat)
	colnames(TCGA.drug.response.mat) = c("TCGA", "Cancer", drugs)
	write.table(TCGA.drug.response.mat, file=paste(ksigmoid,".pred_TCGA.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
	write.table(self.prediction.mat, file=paste(ksigmoid,".pred_CCLE.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
	save(model.list, file=paste(ksigmoid, ".CCLE.model.list.RData", sep=""))
}


