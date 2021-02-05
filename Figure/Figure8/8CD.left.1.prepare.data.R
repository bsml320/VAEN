setwd("/path/to/VAEN/Figure/Figure8")

#################################################################################################################################################################
drug.ccle = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.MIX.pred_TCGA.txt", header=T, as.is=T, sep="\t")

colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, blca.ccle)
}
drug.ccle = cancer.drug.ccle

#################################################################################################################################################################

for(k in 1:length(cancer.types)){
	cancer = cancer.types[k]
	
	if(cancer == "LAML"){
		fn = paste("../../MC3/LAML_wustl", sep="")
	} else {
		fn = paste("../../MC3/",cancer,"_mc3.txt", sep="")
	}
	
	original.maf = read.delim(fn, as.is=T)
	ss = original.maf[,1]; names(ss) = NULL

	tapply(original.maf$gene, original.maf$sample, function(u)unique(as.character(u))) -> ss2gene
	lapply(ss2gene, function(u){sum(u %in% c("BRAF", "KRAS", "NRAS", "HRAS", "EGFR"))}) -> check
	check = unlist(check)
	WT.ss = names(check[which(check==0)])
	match(WT.ss, drug.ccle[,1]) -> ii.2; ii.2 = ii.2[!is.na(ii.2)]
	
	############
	which(original.maf$gene == "BRAF" & grepl("V600", original.maf$Amino_Acid_Change) ) -> ii
	MT.ss = unique(ss[ii])
	match(MT.ss, drug.ccle[,1]) -> ii.1; ii.1 = ii.1[!is.na(ii.1)];
	
	if(cancer == "THCA"){
		subtype = read.delim("../../DATA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
		ss = substr(subtype[,1], 1, 15)
		ss[which(subtype[,2]=="Braf-like")] -> g1
		ss[which(subtype[,2]=="Ras-like")] -> g2

		cbind(drug.ccle[c(ii.1, ii.2), c(1,20,21)], Mut = c(rep("Mut", length(ii.1)), rep("WT", length(ii.2)))) -> tmp
		grp = rep("Other", nrow(tmp))
		grp[which(tmp[,1] %in% g1)] = "Braf"
		grp[which(tmp[,1] %in% g2)] = "Ras"
		tmp = cbind(tmp, grp)
		########################################
		save(tmp, file="THCA.BRAF.plot.RData")
		########################################
	}
	
	if(cancer == "SKCM"){
		cbind(drug.ccle[c(ii.1, ii.2), c(1,20,21)], Mut = c(rep("Mut", length(ii.1)), rep("WT", length(ii.2)))) -> tmp
		########################################
		save(tmp, file="SKCM.BRAF.plot.RData")
		########################################
	}
}
