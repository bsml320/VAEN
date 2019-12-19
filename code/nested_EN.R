### This file was modified based on the Elastic Net method implmented in PrediXcan.
### The full credit should be given to Dr. IM lab: https://github.com/hakyimlab/PrediXcan

generate_fold_ids <- function(n_samples, n_folds=10) {
	n <- ceiling(n_samples / n_folds)
	fold_ids <- rep(1:n_folds, n)
	sample(fold_ids[1:n_samples])
}

calc_R2 <- function(y, y_pred) {
	tss <- sum(y**2)
	rss <- sum((y - y_pred)**2)
	1 - rss/tss
}

calc_corr <- function(y, y_pred) {
	sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

nested_glm = function(x, y, nonzero_idx, n_samples, n_train_test_folds, n_k_folds){
	new.x = as.matrix(x[, nonzero_idx])
	
	R2_folds     <- rep(0, n_train_test_folds)
	corr_folds   <- rep(0, n_train_test_folds)
	zscore_folds <- rep(0, n_train_test_folds)
	pval_folds   <- rep(0, n_train_test_folds)
	
	train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
	for (test_fold in 1:n_train_test_folds) {
		train_idxs <- which(train_test_fold_ids != test_fold)
		test_idxs  <- which(train_test_fold_ids == test_fold)
		x_train    <- new.x[train_idxs, ]
		y_train    <- y[train_idxs]
		x_test     <- new.x[test_idxs, ]
		y_test     <- y[test_idxs]
		
		data = data.frame(y_train, x_train)
		glm(y_train ~ ., data = data) -> fit
		predict(fit, data.frame(x_test)) -> y_pred_test
		
		R2_folds[test_fold] <- calc_R2(y_test, y_pred_test)
		
		corr_folds[test_fold]   <- ifelse(sd(y_pred_test) != 0, cor(y_pred_test, y_test), 0)
		zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
		pval_folds[test_fold]   <- ifelse(sd(y_pred_test) != 0, cor.test(y_pred_test, y_test)$p.value, runif(1))
		cat(test_fold, ".", sep="")
	}
	
	R2_avg <- mean(R2_folds)
	R2_sd  <- sd(R2_folds)
	rho_avg <- mean(corr_folds)
	rho_se  <- sd(corr_folds)
	rho_avg_squared <- rho_avg**2
	
	zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
	zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
	
	# Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
	
	pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
	list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha) {
	# Gets performance estimates for k-fold cross-validated elastic-net models.
	# Splits data into n_train_test_folds disjoint folds, roughly equal in size,
	# and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
	# cross validated. Then get performance measures for how the model predicts on the hold-out
	# fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
	# is there is no correlation between prediction and observed.
	#
	# The mean and standard deviation of R^2 over all folds is then reported, and the p-values
	# are combined using Fisher's method.
	
	R2_folds     <- rep(0, n_train_test_folds)
	corr_folds   <- rep(0, n_train_test_folds)
	zscore_folds <- rep(0, n_train_test_folds)
	pval_folds   <- rep(0, n_train_test_folds)
	# Outer-loop split into training and test set.
	train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
	for (test_fold in 1:n_train_test_folds) {
		train_idxs <- which(train_test_fold_ids != test_fold)
		test_idxs  <- which(train_test_fold_ids == test_fold)
		x_train    <- x[train_idxs, ]
		y_train    <- y[train_idxs]
		x_test     <- x[test_idxs, ]
		y_test     <- y[test_idxs]
		
		# Inner-loop - split up training set for cross-validation to choose lambda.
		cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
		y_pred <- tryCatch({
			  # Fit model with training data.
			  fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids)
			  # Predict test data using model that had minimal mean-squared error in cross validation.
			  predict(fit, x_test, s = 'lambda.min')},
			  # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
			  error = function(cond) rep(mean(y_train), length(y_test)))
		
		R2_folds[test_fold] <- calc_R2(y_test, y_pred)
		# Get p-value for correlation test between predicted y and actual y.
		# If there was no model, y_pred will have var=0, so cor.test will yield NA.
		# In that case, give a random number from uniform distribution, which is what would
		# usually happen under the null.
		corr_folds[test_fold]   <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
		zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
		pval_folds[test_fold]   <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
		cat(test_fold, ".", sep="")
	}
	
	R2_avg <- mean(R2_folds)
	R2_sd  <- sd(R2_folds)
	rho_avg <- mean(corr_folds)
	rho_se  <- sd(corr_folds)
	rho_avg_squared <- rho_avg**2
	
	# Stouffer's method for combining z scores.
	zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
	zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
	
	# Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
	pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
	list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}


main <- function(X, Y, n_folds=10, n_train_test_folds=5, seed=NA, alpha=0.5, null_testing=FALSE, drug) {
	n_samples <- length(Y)
	# Set seed----
	seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
	set.seed(seed)

	model_summary <- c(drug, alpha, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
	weight.mat = c()
	cis_gt = as.matrix(X)
	adj_expression = Y
	
	if (null_testing) adj_expression <- sample(adj_expression)
	
	R2_avg <- NA
	R2_sd <- NA
	pval_est <- NA
	rho_avg <- NA
	rho_se <- NA
	rho_zscore <- NA
	rho_avg_squared <- NA
	zscore_pval <- NA
			
	# Fit on all data
	cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)
	fit <- tryCatch(cv.glmnet(cis_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE),
                    error = function(cond) {message('Error'); message(geterrmessage()); list()}
					)
	
	
	if (length(fit) > 0) {
		cv_R2_folds     <- rep(0, n_folds)
		cv_corr_folds   <- rep(0, n_folds)
		cv_zscore_folds <- rep(0, n_folds)
		cv_pval_folds   <- rep(0, n_folds)
		best_lam_ind    <- which.min(fit$cvm)
		for (j in 1:n_folds) {
			fold_idxs <- which(cv_fold_ids == j)
			adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
			cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
			cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
			cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation
			cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
		}
		
		cv_R2_avg <- mean(cv_R2_folds)
		cv_R2_sd <- sd(cv_R2_folds)
		adj_expr_pred <- predict(fit, as.matrix(cis_gt), s = 'lambda.min')
		training_R2 <- calc_R2(adj_expression, adj_expr_pred)
		
		cv_rho_avg <- mean(cv_corr_folds)
		cv_rho_se <- sd(cv_corr_folds)
		cv_rho_avg_squared <- cv_rho_avg**2
		
		# Stouffer's method for combining z scores.
		cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
		cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
		cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)
        
		beta = coef(fit, fit$lambda.min)
		which(beta!=0) -> nonzero_idx
		nonzero_idx = setdiff(nonzero_idx-1, 0)
		
        if (fit$nzero[best_lam_ind] > 0 & length(nonzero_idx) > 1) {
			set.seed(seed)
			perf_measures <- nested_glm(cis_gt, adj_expression, nonzero_idx, n_samples, n_train_test_folds, n_folds)
			R2_avg <- perf_measures$R2_avg
			R2_sd <- perf_measures$R2_sd
			pval_est <- perf_measures$pval_est
			rho_avg <- perf_measures$rho_avg
			rho_se <- perf_measures$rho_se
			rho_zscore <- perf_measures$rho_zscore
			rho_avg_squared <- perf_measures$rho_avg_squared
			zscore_pval <- perf_measures$zscore_pval
			
			
			weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
			weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
			weight.mat = cbind(weights, weighted_snps)
			
			model_summary <- c(drug, alpha, fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,
                               rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        } else {
          model_summary <- c(drug, alpha, 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,
                             cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                             cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        }
	} else {
		print("fit failed")
		model_summary <- c(drug, alpha, 0, NA, R2_avg, R2_sd, NA, NA, NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, NA, NA, NA, NA, NA, NA)
    }
	
	res.list = list()
	res.list[["model_summary"]] = model_summary
	res.list[["weight.mat"]]    = weight.mat
	res.list[["model"]] = fit
	res.list[["original_Y"]] = Y
	res.list[["self_pred"]] = adj_expr_pred
	res.list
}


