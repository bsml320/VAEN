setwd("/path/to/VAEN/result.EN/dr.GDSC/")

load("GDSC.A.info.RData")

anno = read.delim("../../DATA/GDSC/v17.3_fitted_dose_response.txt", as.is=T)
drugs = sort(unique(anno$DRUG_NAME))

dr.gdsc.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load(paste("01/", best.index, ".GDSC.model.list.RData", sep=""))
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.gdsc.models[[ drug ]] = res.list
}

save(dr.gdsc.models, file="dr.GDSC.A.models.RData")

##########################################################################
##########################################################################
