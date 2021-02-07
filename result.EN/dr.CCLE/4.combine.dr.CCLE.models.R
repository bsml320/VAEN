setwd("/path/to/VAEN/result.EN/dr.CCLE/")

##########################################################################
load("CCLE.A.info.RData")

drugs = colnames(all.F1_R2.mat)
dr.ccle.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("/path/to/VAEN/result.EN/dr.CCLE/01/", best.index ,".CCLE.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.ccle.models[[ drug ]] = res.list
}

save(dr.ccle.models, file="dr.CCLE.A.models.RData")

##########################################################################

load("CCLE.S.info.RData")

dr.ccle.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load( paste("/path/to/VAEN/result.EN/dr.CCLE/01S/", best.index ,".CCLE.model.list.S.RData", sep="") )
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.ccle.models[[ drug ]] = res.list
}

save(dr.ccle.models, file="dr.CCLE.S.models.RData")
