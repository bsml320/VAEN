setwd("/path/to/VAEN/Figure/FigureS10")

###################################################################################################

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

###################################################################################################

BRAF.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
KRAS.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
NRAS.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
HRAS.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
EGFR.mat = matrix(1, nrow=length(drugs), ncol=length(cancer.types))
rownames(BRAF.mat) = drugs; colnames(BRAF.mat) = cancer.types
rownames(KRAS.mat) = drugs; colnames(KRAS.mat) = cancer.types
rownames(NRAS.mat) = drugs; colnames(NRAS.mat) = cancer.types
rownames(HRAS.mat) = drugs; colnames(HRAS.mat) = cancer.types
rownames(EGFR.mat) = drugs; colnames(EGFR.mat) = cancer.types

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
	
	four.ii = c()
	############
	which(original.maf$gene == "BRAF" & original.maf$Amino_Acid_Change == "p.V600E") -> ii
	MT.ss = unique(ss[ii])
	match(MT.ss, drug.ccle[,1]) -> ii.1; ii.1 = ii.1[!is.na(ii.1)]; 
	
	if(length(ii.1) > 10){
		p = c()
		for(k1 in 3:ncol(drug.ccle)){
			BRAF.mat[k1-2, k] = wilcox.test(drug.ccle[ii.1, k1], drug.ccle[ii.2, k1], alternative="greater")$p.value
		}
	}
	
	############
	which(original.maf$gene == "KRAS" & grepl("G12",original.maf$Amino_Acid_Change)) -> ii
	which(original.maf$gene == "KRAS" & grepl("G13",original.maf$Amino_Acid_Change)) -> ii1
	which(original.maf$gene == "KRAS" & grepl("Q61",original.maf$Amino_Acid_Change)) -> ii2
	ii = unique(c(ii, ii1, ii2))
	MT.ss = unique(ss[ii])
	match(MT.ss, drug.ccle[,1]) -> ii.1; ii.1 = ii.1[!is.na(ii.1)]; 
	
	if(length(ii.1) > 10){
		for(k1 in 3:ncol(drug.ccle)){
			KRAS.mat[k1-2, k] = t.test(drug.ccle[ii.1, k1], drug.ccle[ii.2, k1], alternative="greater")$p.value
		}
	}
	
	############
	which(original.maf$gene == "NRAS" & grepl("G12",original.maf$Amino_Acid_Change)) -> ii
	which(original.maf$gene == "NRAS" & grepl("G13",original.maf$Amino_Acid_Change)) -> ii1
	which(original.maf$gene == "NRAS" & grepl("Q61",original.maf$Amino_Acid_Change)) -> ii2
	ii = unique(c(ii, ii1, ii2))
	MT.ss = unique(ss[ii])
	match(MT.ss, drug.ccle[,1]) -> ii.1; ii.1 = ii.1[!is.na(ii.1)]; 
	
	if(length(ii.1) > 10){
		for(k1 in 3:ncol(drug.ccle)){
			NRAS.mat[k1-2, k] = t.test(drug.ccle[ii.1, k1], drug.ccle[ii.2, k1], alternative="greater")$p.value
		}
	}
	
	############
	which(original.maf$gene == "HRAS" & grepl("G12",original.maf$Amino_Acid_Change)) -> ii
	which(original.maf$gene == "HRAS" & grepl("G13",original.maf$Amino_Acid_Change)) -> ii1
	which(original.maf$gene == "HRAS" & grepl("Q61",original.maf$Amino_Acid_Change)) -> ii2
	ii = unique(c(ii, ii1, ii2))
	MT.ss = unique(ss[ii])
	match(MT.ss, drug.ccle[,1]) -> ii.1; ii.1 = ii.1[!is.na(ii.1)]; 
	
	if(length(ii.1) > 10){
		for(k1 in 3:ncol(drug.ccle)){
			HRAS.mat[k1-2, k] = t.test(drug.ccle[ii.1, k1], drug.ccle[ii.2, k1], alternative="greater")$p.value
		}
	}
	
	############
	which(original.maf$gene == "EGFR") -> ii
	grep("746", original.maf$Amino_Acid_Change) -> del.ii
	grep("858", original.maf$Amino_Acid_Change) -> mut.ii
	grep("861", original.maf$Amino_Acid_Change) -> mut2.ii
	grep("719", original.maf$Amino_Acid_Change) -> mut3.ii
	ii = intersect(ii, c(del.ii, mut.ii, mut2.ii, mut3.ii))
	
	MT.ss = unique(ss[ii])
	match(MT.ss, drug.ccle[,1]) -> ii.1; ii.1 = ii.1[!is.na(ii.1)];
	
	if(length(ii.1) > 10){
		for(k1 in 3:ncol(drug.ccle)){
			EGFR.mat[k1-2, k] = t.test(drug.ccle[ii.1, k1], drug.ccle[ii.2, k1], alternative="greater")$p.value
		}
	}
	
}
save(BRAF.mat, EGFR.mat, KRAS.mat, NRAS.mat, HRAS.mat, file="S10.data.RData")

############################################################################################################################
