setwd("/work/RANK.Sigmoid/PCA/dr.GDSC/")

dat = read.table("/work/RANK.Sigmoid/CCLE.4VAE.RANK.tsv", header=T, as.is=T)
res.pca <- prcomp(dat, scale = TRUE)

PPs = res.pca$x[, 1:100]

TCGA = read.table("/work/RANK.Sigmoid/TCGA.4VAE.RANK.tsv", header=T, as.is=T)
ind.sup.coord <- predict(res.pca, newdata = TCGA)

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")

source("/data1_2/jiap/projects/18-CCLE-VAE/new/nested_EN.R")

library("parallel")
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

#####################################################################################
	
	########### Prediction
	original.ss.PP = rownames(PPs)
	sapply(original.ss.PP, function(x){
		new.u = u = strsplit(x, split="\\.")[[1]]
		paste(u[3], u[4], sep="-") -> new.u
		new.u
	}) -> ss.ACH
	names(ss.ACH) = NULL

	
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
		} else {
			res.list = tmp.list[[idx]]
			cat("selected ", idx, ", ", sep="")
		}
		
		
		Ys = matrix(-9, nrow=nrow(PPs), ncol=2)
		rownames(Ys) = rownames(PPs)
		Ys[ c(part1.ach.idx, part2.ach.idx) ,1] = Y
		Ys[ c(part1.ach.idx, part2.ach.idx) ,2] = res.list$self_pred
		res.list$Ys = Ys
		
		model.list[[ drugs[k] ]] = res.list
		
		cat(drugs[k]," end \n", sep="")
	}

	########### match
	save(model.list, file=paste("PCA.GDSC.model.list.RData", sep=""))

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
		
		recall    = cor(Ys[,1], Ys[,2])
		precision = as.numeric(res.list$model_summary[5])
		F1_R2 = 2 * precision * recall/(precision + recall)
		info = rbind(info, c(drug, F1_R2, recall, precision) )
	}
	
colnames(info) = c("Drug", "F1", "Self_PCC", "avg_CV_R2")
write.table(info, file="NOPEER.RANK.PCA.GDSC.txt", row.names=F, col.names=T, quote=F, sep="\t")

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
write.table(TCGA.pred.mat, file="PCA.GDSC.pred_TCGA.txt", quote=F, sep="\t", row.names=FALSE)

self.prediction.mat[,1] = rownames(Ys)
colnames(self.prediction.mat) = c("CELLLINE", drugs)
write.table(self.prediction.mat, file="PCA.GDSC.pred_CCLE.txt", quote=F, sep="\t", row.names=FALSE)

self.pred.mat = cbind(CELLINE=rownames(CCLE.latent.data), CCLE.pred.full.mat)
colnames(self.pred.mat) = c("CELLINE", drugs)
write.table(self.pred.mat, file=paste("PCA.GDSC.pred_CCLE.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)



 
