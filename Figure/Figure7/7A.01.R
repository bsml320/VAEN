#setwd("/path/to/VAEN/Figure/Figure7")

################## prepare cancer mutations

drug.ccle = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.MIX.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

cancer.mutation.list = list()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }

	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	
	if(cancer == "LAML"){
		fn = paste("../../MC3/LAML_wustl", sep="")
	} else {
		fn = paste("../../MC3/",cancer,"_mc3.txt", sep="")
	}
	
	mut.mat = read.delim(fn, as.is=T)
	
	fixed.ss = intersect(blca.ccle[,1], mut.mat[,1] )
	
	mut.mat   = mut.mat[which(mut.mat[,1] %in% fixed.ss),]
	tapply(mut.mat[,1], mut.mat$gene, function(u)setdiff(fixed.ss, u)) -> gene2WT.list
	
	mut.mat = mut.mat[mut.mat$Amino_Acid_Change!="", ]
	mut.mat = mut.mat[mut.mat$effect %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation"), ]
	
	which(mut.mat$effect == "Missense_Mutation") -> ii
	non.mut.mat = mut.mat[-ii, ]
	mis.mut.mat = mut.mat[ii, ]
	
	ii = grepl("deleterious", mis.mut.mat$SIFT) & grepl("damaging", mis.mut.mat$PolyPhen)
	mis.mut.mat = mis.mut.mat[ii, ]
	
	mut.mat = rbind(mis.mut.mat, non.mut.mat)
	
	### with indel
	new.mut.mat.3cols = cbind(gene = mut.mat$gene, sample = mut.mat$sample, AAchange = mut.mat[, "Amino_Acid_Change"])
	new.mut.mat.3cols = new.mut.mat.3cols[new.mut.mat.3cols[,2] %in% fixed.ss, ]
	
	### without indel
	grep("fs", new.mut.mat.3cols[, 3]) -> ii
	noindel.mut.mat.3cols = new.mut.mat.3cols[-ii, ]
	
	mut.list = list()
	mut.list[[ "mut" ]] = new.mut.mat.3cols
	mut.list[[ "mut_wo_indel" ]] = noindel.mut.mat.3cols
	mut.list[[ "gene2WT.list" ]] = gene2WT.list
	mut.list[[ "ss"  ]] = fixed.ss
	cancer.mutation.list[[ cancer ]] = mut.list
}
save(cancer.mutation.list, file="MC3.RData")
print("finish prepare cancer mutation data")

######################################################################

drug.ccle = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.MIX.pred_TCGA.txt", header=T, as.is=T, sep="\t")
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
	new.mut.mat.3cols = mut.list[[ "mut_wo_indel" ]]
	fixed.ss = mut.list[[ "ss"  ]]
	gene2WT.list = mut.list[[ "gene2WT.list" ]]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	
	tapply(new.mut.mat.3cols[,2], new.mut.mat.3cols[,1], function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	print(paste(cancer, length(candidate.genes), nrow(blca.ccle), sep=" "))
	if(length(candidate.genes) < 1)next
	info = rbind(info, c(cancer, length(candidate.genes)))
	
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
save(gene.drug.limma.mat, file=paste("7A.MC3.mut_wo_indel.RData", sep="") )

### LGG
cancer = "LGG"
mut.list = cancer.mutation.list[[ cancer ]]
new.mut.mat.3cols = mut.list[[ "mut_wo_indel" ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]

subtype = read.delim("../../DATA/LGG.subtype.TCGA.txt", as.is=T, header=T)
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
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
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
new.mut.mat.3cols = mut.list[[ "mut_wo_indel" ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]

subtype = read.delim("../../DATA/BRCA.subtype.TCGA.txt", as.is=T, header=T, sep="\t", skip=1)
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
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
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
new.mut.mat.3cols = mut.list[[ "mut_wo_indel" ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]

subtype = read.delim("../../DATA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
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
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
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

###########################################################################

colnames(gene.drug.limma.mat) = c("Cancer", "Gene", "Drug", "p.twosided", "FC", "nMut", "nMutProp")

load("../../DATA/ncbiRefSeq.01252019.hg19.gene2ll.RData")
long.genes = names(gene2ll[gene2ll > 200000])
gene.drug.limma.mat = gene.drug.limma.mat[!gene.drug.limma.mat[,2] %in% long.genes, ]
gene.drug.limma.mat = gene.drug.limma.mat[!gene.drug.limma.mat[,1] %in% c("BRCA", "LGG", "THCA"), ]


gene.drug.limma.mat[which(gene.drug.limma.mat[,3]=="X17.AAG"),3] = "17-AAG"
gsub("\\.", "", gene.drug.limma.mat[,3]) -> new.drug
gene.drug.limma.mat[,3] = new.drug

write.table(gene.drug.limma.mat, file="Figure7A.txt", row.names=F, quote=F, sep="\t")
###########################################################################
gene.drug.limma.mat = read.table("Figure7A.txt", as.is=T, header=T)

###### *****
gene.drug.limma.mat$adjp = NULL
drugs = unique(gene.drug.limma.mat[,3])
cancer.types = unique(gene.drug.limma.mat[,1])
for(k in 1:length(drugs)){
	for(cancer in cancer.types){
		which(gene.drug.limma.mat[,3] == drugs[k] & gene.drug.limma.mat[,1]==cancer) -> ii
		p.adjust(gene.drug.limma.mat[ii,4], method="BH") -> adjp
		gene.drug.limma.mat[ii, "adjp"] = adjp
	}
}
###### *****

which(gene.drug.limma.mat$adjp < 0.005 & gene.drug.limma.mat$FC > 0) -> ii.1
which(gene.drug.limma.mat$adjp < 0.01 & gene.drug.limma.mat$FC < 0) -> ii.2
cc = rep(0, nrow(gene.drug.limma.mat)); cc[ii.1] = 2; cc[ii.2] = 1
gene.drug.limma.mat$color = as.factor(cc)

write.table(gene.drug.limma.mat[which(cc!=0), ], file="Figure7B.txt", row.names=F, quote=F, sep="\t")
write.table(gene.drug.limma.mat, file="Figure7A.txt", row.names=F, quote=F, sep="\t")

