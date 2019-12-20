setwd("/work/result.EN/dr.CCLE/04-mix/01.dSNP/")

################## prepare cancer mutations

cancer.mutation.list = list()
ccle = read.table(file="/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
ccle.ss = gsub("\\.","-",ccle[,1])
ccle[,1] = ccle.ss	

cancer.types = unique(ccle[,2])
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }

	blca.ccle = ccle[which(ccle[,2] == cancer), ]
	sample.type = substr(blca.ccle[,1], 14, 15)
	blca.ccle = blca.ccle[which(sample.type == type.code), ]
	
	fn = paste("/work/data/TCGA/TCGA_annotation/",cancer,"_broad.SNP.avinput.txt.hg19_multianno.txt.hg19_multianno.txt", sep="")
	if(!file.exists(fn)){
		fn = paste("/work/data/TCGA/TCGA_annotation/",cancer,"_bcm.SNP.avinput.txt.hg19_multianno.txt.hg19_multianno.txt", sep="")
		if(cancer == "BRCA" | cancer == "LAML") fn = paste("/work/data/TCGA/TCGA_annotation/",cancer,"_wustl.SNP.avinput.txt.hg19_multianno.txt.hg19_multianno.txt", sep="")
		if(!file.exists(fn)){
			break
		}
	}
	
	header = read.delim(fn, as.is=T, nrows=2, header=F)
	mut.mat = read.delim(fn, as.is=T, skip=2, header=F)
	colnames(mut.mat)[1:126] = header[2, 1:126]
	
	fixed.ss = intersect(blca.ccle[,1], mut.mat[,126] )
	
	mut.mat   = mut.mat[which(mut.mat[,126] %in% fixed.ss),]
	tapply(mut.mat[,126], mut.mat$Gene.refGene, function(u)setdiff(fixed.ss, u)) -> gene2WT.list
	
	mut.mat = mut.mat[which(mut.mat$ExonicFunc.refGene %in% c("nonsynonymous SNV", "stopgain", "stoploss") ), ]
	which(mut.mat$ExonicFunc.refGene == "nonsynonymous SNV") -> idx.1
	which(mut.mat$ExonicFunc.refGene != "nonsynonymous SNV") -> idx.2
	mis.mut.mat = mut.mat[idx.1, ]
	
	grep("pred", colnames(mut.mat)) -> ii
	print(colnames(mut.mat)[ii])
	apply(mis.mut.mat[, ii], 1, function(u)sum(u=="D")) -> check
	del.mut.mat = mis.mut.mat[which(check >= 2), ]
	print(c(cancer, nrow(mis.mut.mat), nrow(del.mut.mat)))
	
	new.mut.mat = rbind(del.mut.mat, mut.mat[idx.2, ])
	new.mut.mat.3cols = cbind(gene = new.mut.mat$Gene.refGene, sample = new.mut.mat[, 126], AAchange = new.mut.mat[, "AAChange.refGene"])
	new.mut.mat.3cols = new.mut.mat.3cols[new.mut.mat.3cols[,2] %in% fixed.ss, ]
	
	mut.list = list()
	mut.list[[ "mut" ]] = new.mut.mat.3cols
	mut.list[[ "gene2WT.list" ]] = gene2WT.list
	mut.list[[ "ss"  ]] = fixed.ss
	cancer.mutation.list[[ cancer ]] = mut.list
	cat(cancer, ".", sep="")
}
save(cancer.mutation.list, file="/work/result.EN/dr.CCLE/04-mix/01.dSNP/cancer.mutation.list.dSNP.WT.RData")
print("finish prepare cancer mutation data")

###################################

######################################################################
library(parallel)
myfun = function(kgene, candidate.genes, new.mut.mat.3cols, gene2WT.list, new.ccle, fixed.ss){
	MT.ss = unique(new.mut.mat.3cols[which(new.mut.mat.3cols[, 1] == candidate.genes[kgene]), 2])
	match(MT.ss, fixed.ss) -> MT.ii
	
	WT.ss = gene2WT.list[[ candidate.genes[kgene] ]]
	match(WT.ss, fixed.ss) -> WT.ii
	
	apply(new.ccle, 2, function(u){
		wilcox.test(u[MT.ii], u[WT.ii])$p.value
	}) -> ps.twosided
	
	apply(new.ccle, 2, function(u){
		mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
	}) -> fcs
	
	res.list = list()
	res.list[["ps.twosided"]] = ps.twosided
	res.list[["fcs"]] = fcs
	res.list[["MT"]] = length(MT.ss)
	return(res.list)
}
######################################################################

drug.ccle = read.table(file="/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

######################################################################
info = c()
gene.drug.limma.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	mut.list = cancer.mutation.list[[ cancer ]]
	new.mut.mat.3cols = mut.list[[ "mut" ]]
	fixed.ss = mut.list[[ "ss"  ]]
	gene2WT.list = mut.list[[ "gene2WT.list" ]]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	
	tapply(new.mut.mat.3cols[,2], new.mut.mat.3cols[,1], function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	print(paste(cancer, length(candidate.genes), nrow(blca.ccle), sep=" "))
	if(length(candidate.genes) < 1)next
	info = rbind(info, c(cancer, length(candidate.genes)))
	
	if(length(candidate.genes) > 500){
		N = ceiling(length(candidate.genes)/30)
		print(paste(N, "steps", sep=" "))
		for(istep in 1:N){
			start = (istep-1) * 30 + 1
			end = istep * 30
			if(end > length(candidate.genes))end = length(candidate.genes)
			mclapply(start:end, myfun, candidate.genes, new.mut.mat.3cols, gene2WT.list, new.ccle, fixed.ss, mc.cores=30) -> test
			for(ks in 1:length(test)){
				kgene = (istep-1) * 30 + ks
				gene.drug.limma.mat = rbind(gene.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, test[[ks]]$ps.twosided, test[[ks]]$fcs, test[[ks]]$MT, test[[ks]]$MT/nrow(blca.ccle)))
			}
			cat(istep, ".", sep="")
		}
	} else {
		for(kgene in 1:length(candidate.genes)){
			ii = which(new.mut.mat.3cols[, 1] == candidate.genes[kgene])
			MT.ss = unique(new.mut.mat.3cols[ii, 2])
			match(MT.ss, fixed.ss) -> MT.ii
			
			WT.ss = gene2WT.list[[ candidate.genes[kgene] ]]
			match(WT.ss, fixed.ss) -> WT.ii
			
			apply(new.ccle, 2, function(u){
				wilcox.test(u[MT.ii], u[WT.ii])$p.value
			}) -> ps.twosided
			
			apply(new.ccle, 2, function(u){
				mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
			}) -> fcs
			
			gene.drug.limma.mat = rbind(gene.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, ps.twosided, fcs, length(MT.ss), length(MT.ss)/nrow(blca.ccle) ))
		}
	}
}
save(gene.drug.limma.mat, file=paste("V15.04.01.dSNP.RData", sep="") )
write.table(info, file="04.01.info.txt", row.names=F, col.names=F, quote=F, sep="\t")

### LGG
cancer = "LGG"
mut.list = cancer.mutation.list[[ cancer ]]
new.mut.mat.3cols = mut.list[[ "mut" ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]

blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == "01"), ]
blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle

subtype = read.delim("/work/data/TCGA/LGG.subtype.TCGA.txt", as.is=T, header=T)
ss = substr(fixed.ss, 1, 12)
subtype[which(subtype[,9]=="coc1"),1] -> g1
subtype[which(subtype[,9]=="coc2"),1] -> g2
subtype[which(subtype[,9]=="coc3"),1] -> g3
subtype.list = list()
subtype.list[["LGG-coc1"]] = g1
subtype.list[["LGG-coc2"]] = g2
subtype.list[["LGG-coc3"]] = g3
	
for(k in 1:length(subtype.list)){
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	substr(new.mut.mat.3cols[,2], 1, 12) -> new.mut.mat.3cols.ss
	subtype.mut.mat = new.mut.mat.3cols[new.mut.mat.3cols.ss %in% subtype.ss, ]
	fixed.ss = unique(subtype.mut.mat[,2])
	
	tapply(subtype.mut.mat[,2], subtype.mut.mat[,1], function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == "LGG" & sample.type == "01"), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	print(paste(cancer, length(candidate.genes), nrow(blca.ccle), sep=" "))
	if(length(candidate.genes) < 1)next
	info = rbind(info, c(cancer, length(candidate.genes)))
	
	for(kgene in 1:length(candidate.genes)){
		ii = which(subtype.mut.mat[, 1] == candidate.genes[kgene])
		MT.ss = unique(subtype.mut.mat[ii, 2])
		match(MT.ss, fixed.ss) -> MT.ii
		
		WT.ss = intersect(gene2WT.list[[ candidate.genes[kgene] ]], fixed.ss)
		match(WT.ss, fixed.ss) -> WT.ii
		
		apply(new.ccle, 2, function(u){
			wilcox.test(u[MT.ii], u[WT.ii])$p.value
		}) -> ps.twosided
		
		apply(new.ccle, 2, function(u){
			mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
		}) -> fcs
		
		gene.drug.limma.mat = rbind(gene.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, ps.twosided, fcs, length(MT.ss), length(MT.ss)/nrow(blca.ccle) ))
	}
}

### BRCA
cancer = "BRCA"
mut.list = cancer.mutation.list[[ cancer ]]
new.mut.mat.3cols = mut.list[[ "mut" ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]

blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == "01"), ]
blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle

subtype = read.delim("/work/data/TCGA/BRCA.subtype.TCGA.txt", as.is=T, header=T, sep="\t", skip=1)
ss = substr(fixed.ss, 1, 12)
subtype[which(subtype$PAM50.mRNA=="Basal-like"),1] -> g1
subtype[which(subtype$PAM50.mRNA=="HER2-enriched"),1] -> g2
subtype[which(subtype$PAM50.mRNA=="Luminal A"),1] -> g3
subtype[which(subtype$PAM50.mRNA=="Luminal B"),1] -> g4
subtype.list = list()
subtype.list[["BRCA-Basal"]] = g1
subtype.list[["BRCA-HER2"]] = g2
subtype.list[["BRCA-LumA"]] = g3
subtype.list[["BRCA-LumB"]] = g4

for(k in 1:length(subtype.list)){
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	substr(new.mut.mat.3cols[,2], 1, 12) -> new.mut.mat.3cols.ss
	subtype.mut.mat = new.mut.mat.3cols[new.mut.mat.3cols.ss %in% subtype.ss, ]
	fixed.ss = unique(subtype.mut.mat[,2])
	
	tapply(subtype.mut.mat[,2], subtype.mut.mat[,1], function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == "BRCA" & sample.type == "01"), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	print(paste(cancer, length(candidate.genes), nrow(blca.ccle), sep=" "))
	if(length(candidate.genes) < 1)next
	info = rbind(info, c(cancer, length(candidate.genes)))
	
	for(kgene in 1:length(candidate.genes)){
		ii = which(subtype.mut.mat[, 1] == candidate.genes[kgene])
		MT.ss = unique(subtype.mut.mat[ii, 2])
		match(MT.ss, fixed.ss) -> MT.ii
		
		WT.ss = intersect(gene2WT.list[[ candidate.genes[kgene] ]], fixed.ss)
		match(WT.ss, fixed.ss) -> WT.ii
		
		apply(new.ccle, 2, function(u){
			wilcox.test(u[MT.ii], u[WT.ii])$p.value
		}) -> ps.twosided
		
		apply(new.ccle, 2, function(u){
			mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
		}) -> fcs
		
		gene.drug.limma.mat = rbind(gene.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, ps.twosided, fcs, length(MT.ss), length(MT.ss)/nrow(blca.ccle) ))
	}
}

### THCA
cancer = "THCA"
mut.list = cancer.mutation.list[[ cancer ]]
new.mut.mat.3cols = mut.list[[ "mut" ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]

blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == "01"), ]
blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle

subtype = read.delim("/work/data/TCGA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
ss = substr(subtype[,1], 1, 15)
ss[which(subtype[,2]=="Braf-like")] -> g1
ss[which(subtype[,2]=="Ras-like")] -> g2
subtype.list = list()
subtype.list[["THCA-Braf-like"]] = g1
subtype.list[["THCA-Ras-like"]] = g2


for(k in 1:length(subtype.list)){
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	subtype.mut.mat = new.mut.mat.3cols[new.mut.mat.3cols[,2] %in% subtype.ss, ]
	fixed.ss = unique(subtype.mut.mat[,2])
	
	tapply(subtype.mut.mat[,2], subtype.mut.mat[,1], function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == "THCA" & sample.type == "01"), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	print(paste(cancer, length(candidate.genes), nrow(blca.ccle), sep=" "))
	if(length(candidate.genes) < 1)next
	info = rbind(info, c(cancer, length(candidate.genes)))
	
	for(kgene in 1:length(candidate.genes)){
		ii = which(subtype.mut.mat[, 1] == candidate.genes[kgene])
		MT.ss = unique(subtype.mut.mat[ii, 2])
		match(MT.ss, fixed.ss) -> MT.ii
		
		WT.ss = intersect(gene2WT.list[[ candidate.genes[kgene] ]], fixed.ss)
		match(WT.ss, fixed.ss) -> WT.ii
		
		apply(new.ccle, 2, function(u){
			wilcox.test(u[MT.ii], u[WT.ii])$p.value
		}) -> ps.twosided
		
		apply(new.ccle, 2, function(u){
			mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
		}) -> fcs
		
		gene.drug.limma.mat = rbind(gene.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, ps.twosided, fcs, length(MT.ss), length(MT.ss)/nrow(blca.ccle) ))
	}
}

save(gene.drug.limma.mat, file=paste("04.01.dSNP.RData", sep="") )
write.table(info, file="04.01.info.txt", row.names=F, col.names=F, quote=F, sep="\t")


