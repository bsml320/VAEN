#setwd("/path/to/VAEN/Figure/Figure3")

#####################################################################################

load("../../result.EN/dr.GDSC/01/1.GDSC.model.list.RData")
gdsc.model.list = model.list
gdsc.obsd = c()
for(k in 1:length(gdsc.model.list)){
	res.list = gdsc.model.list[[k]]
	drug = names(gdsc.model.list)[k]
	Ys = res.list$Ys
	
	gdsc.obsd = cbind(gdsc.obsd, Ys[,1])
}
colnames(gdsc.obsd) = names(gdsc.model.list)

#######################################

gdsc.pred      = read.table(file="../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_GDSC.txt", header=T, as.is=T, sep="\t")
gdsc.pred.full = read.table(file="../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_GDSC.full.txt", header=T, as.is=T, sep="\t")

colnames(gdsc.pred)      = c("Sample", colnames(gdsc.obsd))
colnames(gdsc.pred.full) = c("Sample", colnames(gdsc.obsd))

#####################################################################################
########### Prediction
original.ss.PP = rownames(gdsc.obsd)
sapply(original.ss.PP, function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
names(tt) = NULL

#####################################################################################

anno = read.delim("../../DATA/GDSC/v17.3_fitted_dose_response.txt", as.is=T)
drugs = sort(unique(anno$DRUG_NAME))

#####################################################################################

full.drug.by.tissue.mat = pred.drug.by.tissue.mat = obsd.drug.by.tissue.mat = matrix(1, nrow=length(drugs), ncol=length(unique(tt)))
rownames(obsd.drug.by.tissue.mat) = rownames(pred.drug.by.tissue.mat) = rownames(full.drug.by.tissue.mat) = drugs
colnames(obsd.drug.by.tissue.mat) = colnames(pred.drug.by.tissue.mat) = colnames(full.drug.by.tissue.mat) = sort(unique(tt))

for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	Ys = cbind(gdsc.obsd[, drug], gdsc.pred[, drug], gdsc.pred.full[, drug])
	
	sapply(rownames(Ys), function(x){
		strsplit(x, split="\\.")[[1]][1] -> u
		strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
	}) -> tt
	names(tt) = NULL
	
	for(kt in 1:ncol(obsd.drug.by.tissue.mat)){
		x.tissue = colnames(obsd.drug.by.tissue.mat)[kt]
		
		X = ifelse(tt %in% x.tissue, 1, 0)
		
		which(Ys[,1]!=-9) -> ii
		Y1 = Ys[ii, 1]
		X1 = X[ii]
		sum(X1==1) -> check
		label = ifelse(mean(Y1[which(X1==1)]) > mean(Y1[which(X1==0)]), 1, -1)
		if(check > 5)obsd.drug.by.tissue.mat[kdrug, kt] = label * t.test(Y1[which(X1==1)], Y1[which(X1==0)])$p.value
		
		which(Ys[,2]!=-9) -> ii
		Y2 = Ys[ii, 2]
		X2 = X[ii]
		sum(X2==1) -> check
		label = ifelse(mean(Y2[which(X2==1)]) > mean(Y2[which(X2==0)]), 1, -1)
		if(check > 5)pred.drug.by.tissue.mat[kdrug, kt] = label * t.test(Y2[which(X2==1)], Y2[which(X2==0)])$p.value
		
		Y3 = Ys[,3]
		label = ifelse(mean(Y3[which(X==1)]) > mean(Y3[which(X==0)]), 1, -1)
		full.drug.by.tissue.mat[kdrug, kt] = label * t.test(Y3[which(X==1)], Y3[which(X==0)])$p.value
	}
}


apply(obsd.drug.by.tissue.mat, 2, var) -> obsd.check
apply(pred.drug.by.tissue.mat, 2, var) -> pred.check

print(which(obsd.check==0))
print(which(pred.check==0))

obsd.drug.by.tissue.mat = obsd.drug.by.tissue.mat[, which(obsd.check != 0)]
pred.drug.by.tissue.mat = pred.drug.by.tissue.mat[, which(obsd.check != 0)]
full.drug.by.tissue.mat = full.drug.by.tissue.mat[, which(obsd.check != 0)]

########################################################################################################################
predicted.dr = read.table(paste("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", sep=""), header=T, as.is=T, sep="\t")

drugs = colnames(predicted.dr)[c(-1, -2)]
cancer.types = unique(predicted.dr[,2])
sample.type = substr(predicted.dr[,1], 14, 15)

cancer.predicted.dr = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = predicted.dr[which(predicted.dr[,2] == cancer & sample.type == type.code), ]
	cancer.predicted.dr = rbind(cancer.predicted.dr, blca.ccle)
}

TCGA.drug.by.cancer.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(TCGA.drug.by.cancer.mat) = drugs
colnames(TCGA.drug.by.cancer.mat) = sort(cancer.types)

for(kdrug in 1:length(drugs)){
	Y = cancer.predicted.dr[, drugs[kdrug]]
	
	for(kt in 1:ncol(TCGA.drug.by.cancer.mat)){
		x.cancer = colnames(TCGA.drug.by.cancer.mat)[kt]
		X = ifelse(cancer.predicted.dr[,2] %in% x.cancer, 1, 0)
		if(sum(X==1) < 5)next
		label = -1
		if(mean(Y[which(X==1)]) > mean(Y[which(X==0)])) label = 1
		
		p = t.test(Y[which(X==1)], Y[which(X==0)])$p.value
		if( p < 1e-100) p = 1e-101
		
		TCGA.drug.by.cancer.mat[kdrug, kt] = label * p
	}
}

########################################################################################################################
########################################################################################################################
predicted.dr = read.table(paste("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", sep=""), header=T, as.is=T, sep="\t")

immune.cancer = c("LAML", "DLBC", "THYM")
which(predicted.dr[,2] %in% immune.cancer) -> ii
predicted.dr = predicted.dr[-ii, ]

drugs = colnames(predicted.dr)[c(-1, -2)]
cancer.types = unique(predicted.dr[,2])
sample.type = substr(predicted.dr[,1], 14, 15)

cancer.predicted.dr = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = predicted.dr[which(predicted.dr[,2] == cancer & sample.type == type.code), ]
	cancer.predicted.dr = rbind(cancer.predicted.dr, blca.ccle)
}

TCGA.sensitive.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(TCGA.sensitive.mat) = drugs
colnames(TCGA.sensitive.mat) = sort(cancer.types)

sensitive.prop = c()
for(kdrug in 1:length(drugs)){
	Y = cancer.predicted.dr[, drugs[kdrug]]
	Y > quantile(Y, probs=.95) -> sensitive.ii
	for(kt in 1:ncol(TCGA.sensitive.mat)){
		x.cancer = colnames(TCGA.sensitive.mat)[kt]
		X = cancer.predicted.dr[,2] %in% x.cancer
		if(kdrug == which(drugs == "PLX.4720"))sensitive.prop = c(sensitive.prop, sum(X & sensitive.ii)/sum(X) )
		
		if(sum(X & sensitive.ii) < 5)next
		
		table(sensitive.ii, X) -> mm
		mm[2:1, 2:1] -> new.mm
		TCGA.sensitive.mat[kdrug, kt] = fisher.test(new.mm)$p.value
	}
}
names(sensitive.prop) = colnames(TCGA.sensitive.mat)

########################################################################################################################
TCGA.resistant.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(TCGA.resistant.mat) = drugs
colnames(TCGA.resistant.mat) = sort(cancer.types)
resistant.prop = c()

for(kdrug in 1:length(drugs)){
	Y = cancer.predicted.dr[, drugs[kdrug]]
	Y < quantile(Y, probs=.05) -> resistant.ii
	for(kt in 1:ncol(TCGA.sensitive.mat)){
		x.cancer = colnames(TCGA.sensitive.mat)[kt]
		X = cancer.predicted.dr[,2] %in% x.cancer
		if(kdrug == which(drugs == "PLX.4720"))resistant.prop = c(resistant.prop, sum(X & resistant.ii)/sum(X) )
		if(sum(X & resistant.ii) < 5)next
		
		table(resistant.ii, X) -> mm
		mm[2:1, 2:1] -> new.mm
		TCGA.resistant.mat[kdrug, kt] = fisher.test(new.mm)$p.value
	}
}
names(resistant.prop) = colnames(TCGA.sensitive.mat)

save(resistant.prop, sensitive.prop, obsd.drug.by.tissue.mat, pred.drug.by.tissue.mat, full.drug.by.tissue.mat, TCGA.drug.by.cancer.mat, TCGA.sensitive.mat, TCGA.resistant.mat, file="3EF.data.RData")

write.table(TCGA.sensitive.mat, file="Figure3F.sensitive.txt", quote=F)
write.table(TCGA.resistant.mat, file="Figure3F.resistant.txt", quote=F)
