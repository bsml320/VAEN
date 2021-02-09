setwd("/path/to/VAEN/Figure/FigureS4/")

anno = read.csv("../../DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))

calc_R2 <- function(y, y_pred) {
	tss <- sum(y**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

calc_R3 <- function(y, y_pred) {
	tss <- sum( (y-mean(y))**2 )
	rss <- sum( (y - y_pred)**2)
	1 - rss/tss
}

calc_MSE <- function(y, y_pred) {
	mean( (y - y_pred)^2 )
}

generate_fold_ids <- function(n_samples, n_folds=10) {
	n <- ceiling(n_samples / n_folds)
	fold_ids <- rep(1:n_folds, n)
	sample(fold_ids[1:n_samples])
}

#####################################################################################

all.ind_CV_R3.mat = all.in_sample_R2.mat = all.avg_CV_R2.mat = all.avg_CV_R3.mat = all.in_sample_MSE.mat = matrix(0, nrow=100, ncol=length(drugs))
colnames(all.ind_CV_R3.mat) = colnames(all.in_sample_R2.mat) = colnames(all.avg_CV_R2.mat) = colnames(all.in_sample_MSE.mat) = colnames(all.avg_CV_R3.mat) = drugs

for(ksigmoid in 1:100){
	load(paste("../../result.EN/dr.CCLE/01/", ksigmoid, ".CCLE.model.list.RData", sep=""))
	
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		model.list[[ drugs[kdrug] ]] -> res.list
		if(length(res.list) == 0)next
		fit <- res.list$model
		Ys = res.list$Ys
		which(Ys[,1]!=-9) -> ii
		Ys = Ys[ii, ]
		if(sd(Ys[,2]) == 0)next
		
		all.in_sample_R2.mat[ksigmoid, kdrug] = cor(Ys[,1], Ys[,2])
		all.avg_CV_R2.mat[ksigmoid, kdrug]    = as.numeric(res.list$model_summary[5])
		all.in_sample_MSE.mat[ksigmoid, kdrug] = as.numeric(res.list$model_summary[7])
		
		cv_fold_ids = fit$foldid
		best_lam_ind    <- which.min(fit$cvm)
		cv_R2_folds = rep(0, 10)
		for(j in 1:10){
			fold_idxs <- which(cv_fold_ids == j)
			adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
			cv_R2_folds[j] <- calc_R3(res.list$original_Y[fold_idxs], adj_expr_fold_pred)
		}
		all.avg_CV_R3.mat[ksigmoid, kdrug]    = mean(cv_R2_folds)
		
		all.ind_CV_R3.mat[ksigmoid, kdrug] = as.numeric(res.list$model_summary[5]) * as.numeric(res.list$model_summary[7])/(as.numeric(res.list$model_summary[5])+as.numeric(res.list$model_summary[7]))
	}
	cat(ksigmoid, ".", sep="")
}

#################################################################################################################################

info = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	tmp = cbind(idx = c(1:100), all.in_sample_MSE.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug], all.avg_CV_R3.mat[, drug], all.ind_CV_R3.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),] ### avg_CV_R2
	
	info = rbind(info, c(drug, tmp[1,]) )
}
colnames(info) = c("Drug", "Run", "Self_MSE", "Self_PCC", "avg_CV_R2", "avg_CV_R3", "ind_CV_R3")

write.table(info, file="VAE.best.info.txt", row.names=F, col.names=T, quote=F, sep="\t")
save(all.in_sample_R2.mat, all.avg_CV_R2.mat, all.avg_CV_R3.mat, file="info.RData")
