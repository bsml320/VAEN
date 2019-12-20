setwd("/work/result.EN/dr.GDSC/04-mix/12.TML")

##################################################################################################
drug.ccle = read.table("/work/result.EN/dr.GDSC/01/GDSC.MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
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

###################################################################################################################

cancer.drug.tml.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer), ]
	
	fn = paste("/work/data/TCGA/",cancer, "/mutation_broad", sep="")
	if(!file.exists(fn)){
		fn = paste("/work/data/TCGA/",cancer, "/mutation_bcm", sep="")
		if(cancer == "BRCA" | cancer == "LAML") fn = paste("/work/data/TCGA/",cancer, "/mutation_wustl", sep="")
		if(!file.exists(fn)){
			break
		}
	}
	mut.mat = read.delim(fn, header=T, as.is=T)
	
	fixed.ss = intersect(blca.ccle[,1], mut.mat[,1] )
	if(length(fixed.ss) < 50){
		print(c("skip ", cancer))
		next
	}
	
	tapply(mut.mat[,1], mut.mat[,1], length) -> sample2gene.length
	sample2gene.length = sample2gene.length[fixed.ss]
	sample2gene.length[which(is.na(sample2gene.length))] = 0
	names(sample2gene.length) = fixed.ss
	
	
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	
	
	for(kdrug in 1:ncol(new.ccle)){
		which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9)) -> ii
		wilcox.test(sample2gene.length[ii], sample2gene.length[-ii])$p.value -> p
		mean(sample2gene.length[ii])/mean(sample2gene.length[-ii]) -> FC
		cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind(cancer, drug = colnames(new.ccle)[kdrug], p=p, FC=FC) )
	}
}
write.table(cancer.drug.tml.mat, file="GDSC.cancer.drug.tml.wilcox.txt", row.names=F, quote=F, sep="\t")

###################################################################################################################

cancer.drug.tml.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer), ]
	
	fn = paste("/work/data/TCGA/",cancer, "/mutation_broad", sep="")
	if(!file.exists(fn)){
		fn = paste("/work/data/TCGA/",cancer, "/mutation_bcm", sep="")
		if(cancer == "BRCA" | cancer == "LAML") fn = paste("/work/data/TCGA/",cancer, "/mutation_wustl", sep="")
		if(!file.exists(fn)){
			break
		}
	}
	mut.mat = read.delim(fn, header=T, as.is=T)
	
	fixed.ss = intersect(blca.ccle[,1], mut.mat[,1] )
	if(length(fixed.ss) < 50){
		print(c("skip ", cancer))
		next
	}
	
	tapply(mut.mat[,1], mut.mat[,1], length) -> sample2gene.length
	sample2gene.length = sample2gene.length[fixed.ss]
	sample2gene.length[which(is.na(sample2gene.length))] = 0
	names(sample2gene.length) = fixed.ss
	
	
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	
	
	for(kdrug in 1:ncol(new.ccle)){
		which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9)) -> ii
		t.test( log(sample2gene.length[ii]+1), log(sample2gene.length[-ii]+1) )$p.value -> p
		mean(sample2gene.length[ii])/mean(sample2gene.length[-ii]) -> FC
		cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind(cancer, drug = colnames(new.ccle)[kdrug], p=p, FC=FC) )
	}
}
write.table(cancer.drug.tml.mat, file="GDSC.cancer.drug.tml.ttest.txt", row.names=F, quote=F, sep="\t")
