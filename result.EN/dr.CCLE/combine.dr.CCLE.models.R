setwd("/path/to/GitHub/result.EN/dr.CCLE/")
load("/path/to/result.EN/dr.CCLE/CCLE.A.info.RData")

anno = read.csv("/path/to/DATA/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))

dr.ccle.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("/path/to/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.ccle.models[[ drug ]] = res.list
}

save(dr.ccle.models, file="dr.CCLE.A.models.RData")


##########################################################################
##########################################################################

load("/path/to/result.EN/dr.CCLE/CCLE.S.info.RData")

dr.ccle.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("/path/to/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.S.RData", sep="") )
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.ccle.models[[ drug ]] = res.list
}

save(dr.ccle.models, file="dr.CCLE.S.models.RData")
