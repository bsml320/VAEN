setwd("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/EXPR/Figure2")
load("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/01/EXPR/EXPR.CCLE.model.list.RData")

anno = read.csv("/data/jiap/data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))

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
		
		info = rbind(info, c(drug, calc_MSE(Ys[,1], Ys[,2]), cor(Ys[,1], Ys[,2]), as.numeric(res.list$model_summary[5]), mean(cv_R2_folds) ) )
	}
	
colnames(info) = c("Drug", "Self_MSE", "Self_PCC", "avg_CV_R2", "avg_CV_R3")

write.table(info, file="EXPR.CCLE.txt", row.names=F, col.names=T, quote=F, sep="\t")
