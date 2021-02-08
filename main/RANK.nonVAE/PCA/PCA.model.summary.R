setwd("/path/to/VAEN/main/RANK.nonVAE/PCA/dr.CCLE/01")

load("PCA.CCLE.model.list.07292020.RData")

info = c()
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		model.list[[ drugs[kdrug] ]] -> res.list
		if(length(res.list) == 0)next
		Ys = res.list$Ys
		which(Ys[,1]!=-9) -> ii
		Ys = Ys[ii, ]
		if(sd(Ys[,2]) == 0)next
		
		fit <- res.list$model
		cv_fold_ids = fit$foldid
		best_lam_ind    <- which.min(fit$cvm)
		cv_R2_folds = rep(0, 10)
		for(j in 1:10){
			fold_idxs <- which(cv_fold_ids == j)
			adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
			cv_R2_folds[j] <- calc_R3(res.list$original_Y[fold_idxs], adj_expr_fold_pred)
		}
		
		info = rbind(info, c(drug, calc_MSE(Ys[,1], Ys[,2]), cor(Ys[,1], Ys[,2]), as.numeric(res.list$model_summary[5]), mean(cv_R2_folds), mean(cv_R2_folds) ) )
	}
	
colnames(info) = c("Drug", "Self_MSE", "Self_PCC", "avg_CV_R2", "avg_CV_R3", "ind_CV_R3")

cbind(vae.info[,"avg_CV_R3"], info[,5])
sum(vae.info[,"avg_CV_R3"] > info[,5])

write.table(info, file="PCA.CCLE.07292020.txt", row.names=F, col.names=T, quote=F, sep="\t")


################################################################################################

load("/path/to/VAEN/main/RANK.nonVAE/PCA/dr.GDSC/01/PCA.GDSC.model.list.07292020.RData")
drugs = names(model.list)
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
		
		cv_fold_ids = fit$foldid
		best_lam_ind    <- which.min(fit$cvm)
		cv_R2_folds = rep(0, 10)
		for(j in 1:10){
			fold_idxs <- which(cv_fold_ids == j)
			adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
			cv_R2_folds[j] <- calc_R3(res.list$original_Y[fold_idxs], adj_expr_fold_pred)
		}
		
		cv_fold_ids = generate_fold_ids(nrow(Ys), 10)
		best_lam_ind    <- which.min(fit$cvm)
		cv_R2_folds = rep(0, 10)
		for(j in 1:10){
			fold_idxs <- which(cv_fold_ids == j)
			adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
			cv_R2_folds[j] <- calc_R3(res.list$original_Y[fold_idxs], adj_expr_fold_pred)
		}
		
		info = rbind(info, c(drug, calc_MSE(Ys[,1], Ys[,2]), cor(Ys[,1], Ys[,2]), as.numeric(res.list$model_summary[5]), mean(cv_R2_folds), mean(cv_R2_folds) ) )
	}
	
colnames(info) = c("Drug", "Self_MSE", "Self_PCC", "avg_CV_R2", "avg_CV_R3", "ind_CV_R3")

write.table(info, file="PCA.GDSC.07292020.txt", row.names=F, col.names=T, quote=F, sep="\t")
