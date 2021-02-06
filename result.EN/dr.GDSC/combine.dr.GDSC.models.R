#setwd("D:/UTH/work/18-VAE/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/GitHub/result.EN/dr.GDSC/")
#setwd("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.GDSC/01/F1-W5-PCC/")

load("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.GDSC/01/F1-W5-PCC/F1-W5-PCC.info.RData")

anno = read.delim("/data1_2/jiap/projects/18-CCLE-VAE/DATA/GDSC//v17.3_fitted_dose_response.txt", as.is=T)
drugs = sort(unique(anno$DRUG_NAME))

dr.gdsc.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	tmp = cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load(paste("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.GDSC/01/", best.index, ".GDSC.model.list.RData", sep=""))
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.gdsc.models[[ drug ]] = res.list
}

save(dr.gdsc.models, file="dr.GDSC.A.models.RData")


##########################################################################
##########################################################################

load("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.GDSC/01S/F1-W5-PCC/S.F1-W5-PCC.info.RData")

dr.gdsc.models = list()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	if(drug == "X17.AAG")drug = "17-AAG"
	gsub("\\.", "-", drug) -> drug
	tmp = cbind(idx = c(1:100), solid.F1_R2.mat[, drug], solid.in_sample_R2.mat[, drug], solid.avg_CV_R2.mat[, drug] )
	tmp = tmp[order(tmp[,4], decreasing=T),]
	best.index = tmp[1,1]
	
	load(paste("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.GDSC/01S/", best.index, ".GDSC.model.list.S.RData", sep=""))
	model.list[[ drug ]] -> res.list
	res.list[[ "best_index" ]] = best.index
	dr.gdsc.models[[ drug ]] = res.list
}

save(dr.gdsc.models, file="dr.GDSC.S.models.RData")
