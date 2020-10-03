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


library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")

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

anno = read.delim("/work/data/GDSC/v17.3_fitted_dose_response.txt", as.is=T)
drugs = sort(unique(anno$DRUG_NAME))
cell.line.anno = read.csv("/work/data/DepMap-2018q3-celllines.csv", as.is=T)

for(ksigmoid in start:end){
	
	cat("ksigmoid = ", ksigmoid, "\n", sep="")

	#############
	model_summary_file <- paste(ksigmoid, ".model_summary.txt", sep="")
	model_summary_cols <- c('Drug', 'alpha', 'n_snps_in_model', 'lambda_min_mse',
							  'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
							  'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
							  'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
	
	write(model_summary_cols, file = model_summary_file, ncol = 21, sep = '\t')
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

	original.ss.PP = rownames(PPs)
	sapply(original.ss.PP, function(x){
		new.u = u = strsplit(x, split="\\.")[[1]]
		paste(u[3], u[4], sep="-") -> new.u
		new.u
	}) -> ss.ACH
	names(ss.ACH) = NULL

	########### original drug data
	self.prediction.mat = matrix(-9, nrow=length(unique(anno[,4])), ncol=length(drugs)+2)
	self.prediction.mat[,1] = unique(anno[,"COSMIC_ID"])
	self.prediction.mat[,2] = unique(anno[,"CELL_LINE_NAME"])
	colnames(self.prediction.mat) = c("CELLLINE", "Type", drugs)
	

	TCGA.drug.response.mat = c()
	model.list = list()
	for(k in 1:length(drugs)){
		cat(drugs[k], " ======== start\n", sep="")
		anno.1 = anno[which(anno$DRUG_NAME==drugs[k]),]
	
		match(anno.1$COSMIC_ID, cell.line.anno$COSMIC_ID) -> idx
		
		anno.1.match1 = anno.1[which(  anno.1$COSMIC_ID %in% cell.line.anno$COSMIC_ID ), ]
		gdsc.anno.1.match2 = anno.1[which(  !anno.1$COSMIC_ID %in% cell.line.anno$COSMIC_ID), ]
		
		match(anno.1.match1$COSMIC_ID, cell.line.anno$COSMIC_ID) -> idx1
		match(cell.line.anno[idx1, 1], ss.ACH) -> part1.ach.idx
		anno.1.match1 = anno.1.match1[which(!is.na(part1.ach.idx)),]
		match(anno.1.match1$COSMIC_ID, cell.line.anno$COSMIC_ID) -> idx1
		match(cell.line.anno[idx1, 1], ss.ACH) -> part1.ach.idx
		Y = -anno.1.match1[, "LN_IC50"]
		
		sapply(gdsc.anno.1.match2$CELL_LINE_NAME, function(u){
			gsub("-", "", u)
		}) -> gdsc.cell.name
		union(which(gdsc.cell.name %in% cell.line.anno$Aliases), which(gdsc.anno.1.match2$CELL_LINE_NAME %in% cell.line.anno$Aliases)) -> ii
		gdsc.anno.1.match2 = gdsc.anno.1.match2[ii, ]
		gdsc.cell.name = gdsc.cell.name[ii]
		
		idx2 = c()
		for(kk in 1:nrow(gdsc.anno.1.match2)){
			if(gdsc.anno.1.match2$CELL_LINE_NAME[kk] %in% cell.line.anno$Aliases) {
				idx2 = c(idx2, match(gdsc.anno.1.match2$CELL_LINE_NAME[kk], cell.line.anno$Aliases))
			} else {
				idx2 = c(idx2, match(gdsc.cell.name[kk], cell.line.anno$Aliases))
			}
		}
		
		match(cell.line.anno[idx2, 1], ss.ACH) -> part2.ach.idx
		gdsc.anno.1.match2 = gdsc.anno.1.match2[which(!is.na(part2.ach.idx)),]
		
		idx2 = c()
		for(kk in 1:nrow(gdsc.anno.1.match2)){
			if(gdsc.anno.1.match2$CELL_LINE_NAME[kk] %in% cell.line.anno$Aliases) {
				idx2 = c(idx2, match(gdsc.anno.1.match2$CELL_LINE_NAME[kk], cell.line.anno$Aliases))
			} else {
				idx2 = c(idx2, match(gsub("-", "", gdsc.anno.1.match2$CELL_LINE_NAME[kk]), cell.line.anno$Aliases))
			}
		}
		
		match(cell.line.anno[idx2, 1], ss.ACH) -> part2.ach.idx
		new.Y = -gdsc.anno.1.match2[, "LN_IC50"]
		names(new.Y) = gdsc.anno.1.match2[, "CELL_LINE_NAME"]
		Y = c(Y, new.Y)
		
		train.data = PPs[c(part1.ach.idx, part2.ach.idx), ]
		
		tmp.list = list()
		drug = drugs[k]
		if(drug == "VNLG/124")drug = "VNLG.124"
		
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
		Ys[ c(part1.ach.idx, part2.ach.idx) ,1] = Y
		Ys[ c(part1.ach.idx, part2.ach.idx) ,2] = res.list$self_pred
		res.list$Ys = Ys
		
		model.list[[ drugs[k] ]] = res.list
		
		match(anno.1.match1$COSMIC_ID,self.prediction.mat[,1]) -> pii.1
		match(gdsc.anno.1.match2$CELL_LINE_NAME,self.prediction.mat[,2]) -> pii.2
		self.prediction.mat[c(pii.1, pii.2), k+2] = res.list$self_pred
		
		cat(drugs[k]," end \n", sep="")
	}

	########### match
	gsub("\\.", "-", TCGA.pred[,1]) -> ss; TCGA.pred[,1] = ss
	match(TCGA.pred[,1], TCGA.ss.mat[,1]) -> ii
	TCGA.drug.response.mat = cbind(TCGA.pred[,1], TCGA.ss.mat[ii, 2], TCGA.drug.response.mat)
	colnames(TCGA.drug.response.mat) = c("TCGA", "Cancer", drugs)
	write.table(TCGA.drug.response.mat, file=paste(ksigmoid,".GDSC.pred_TCGA.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
	write.table(self.prediction.mat, file=paste(ksigmoid,".pred_GDSC.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
	save(model.list, file=paste(ksigmoid, ".GDSC.model.list.RData", sep=""))
}


