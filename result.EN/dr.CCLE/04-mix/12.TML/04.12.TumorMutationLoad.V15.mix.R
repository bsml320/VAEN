setwd("/work/result.EN/dr.CCLE/04-mix/12.TML")

##################################################################################################

drug.ccle = read.table(file="/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
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

##################################################################################################

cancer.tml.list = list()
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
	if(length(fixed.ss) < 50){	next }
	
	tapply(mut.mat[,1], mut.mat[,1], length) -> sample2gene.length
	sample2gene.length = sample2gene.length[fixed.ss]
	sample2gene.length[which(is.na(sample2gene.length))] = 0
	names(sample2gene.length) = fixed.ss
	
	cancer.tml.list[[cancer]] = sample2gene.length
}
save(cancer.tml.list, file="cancer.tml.list.RData")

###################################################################################################################

cancer.drug.tml.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
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
	
	if(cancer == "BRCA"){
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
			subtype.ss = subtype.list[[k]]
			substr( names(sample2gene.length), 1, 12) -> old.ss
			fixed.ss = intersect(subtype.ss, old.ss)
			
			blca.ccle = drug.ccle[match(fixed.ss, substr(drug.ccle[,1], 1, 12) ),]
			apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
			new.sample2gene.length = sample2gene.length[match(fixed.ss, old.ss)]
			names(new.sample2gene.length) = fixed.ss
			
			for(kdrug in 1:ncol(new.ccle)){
				which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9)) -> ii
				wilcox.test(new.sample2gene.length[ii], new.sample2gene.length[-ii])$p.value -> p
				mean(new.sample2gene.length[ii])/mean(new.sample2gene.length[-ii]) -> FC
				cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind( names(subtype.list)[k], drug = colnames(new.ccle)[kdrug], p=p, FC=FC) )
			}
		}
	}
	
	if(cancer == "LGG"){
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
			subtype.ss = subtype.list[[k]]
			substr( names(sample2gene.length), 1, 12) -> old.ss
			fixed.ss = intersect(subtype.ss, old.ss)
			
			blca.ccle = drug.ccle[match(fixed.ss, substr(drug.ccle[,1], 1, 12) ),]
			apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
			new.sample2gene.length = sample2gene.length[match(fixed.ss, old.ss)]
			names(new.sample2gene.length) = fixed.ss
			
			for(kdrug in 1:ncol(new.ccle)){
				which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9)) -> ii
				wilcox.test(new.sample2gene.length[ii], new.sample2gene.length[-ii])$p.value -> p
				mean(new.sample2gene.length[ii])/mean(new.sample2gene.length[-ii]) -> FC
				cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind( names(subtype.list)[k], drug = colnames(new.ccle)[kdrug], p=p, FC=FC) )
			}
		}
	}
	
	if(cancer == "THCA"){
		subtype = read.delim("/work/data/TCGA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
		ss = substr(subtype[,1], 1, 15)
		ss[which(subtype[,2]=="Braf-like")] -> g1
		ss[which(subtype[,2]=="Ras-like")] -> g2
		subtype.list = list()
		subtype.list[["THCA-Braf-like"]] = g1
		subtype.list[["THCA-Ras-like"]] = g2
		
		for(k in 1:length(subtype.list)){
			subtype.ss = subtype.list[[k]]
			fixed.ss = intersect(subtype.ss, names(sample2gene.length))
			
			blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1] ),]
			apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
			new.sample2gene.length = sample2gene.length[fixed.ss]
			
			for(kdrug in 1:ncol(new.ccle)){
				which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9)) -> ii
				wilcox.test(new.sample2gene.length[ii], new.sample2gene.length[-ii])$p.value -> p
				mean(new.sample2gene.length[ii])/mean(new.sample2gene.length[-ii]) -> FC
				cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind( names(subtype.list)[k], drug = colnames(new.ccle)[kdrug], p=p, FC=FC) )
			}
		}
	}
}

write.table(cancer.drug.tml.mat, file="CCLE.cancer.drug.tml.wilcox.txt", row.names=F, quote=F, sep="\t")

###################################################################################################################

cancer.drug.tml.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
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
		t.test(sample2gene.length[ii], sample2gene.length[-ii])$p.value -> p
		mean(sample2gene.length[ii])/mean(sample2gene.length[-ii]) -> FC
		cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind(cancer, drug = colnames(new.ccle)[kdrug], p=p, FC=FC) )
	}
}
write.table(cancer.drug.tml.mat, file="CCLE.cancer.drug.tml.ttest.txt", row.names=F, quote=F, sep="\t")

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
		new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9) -> cat.a
		sample2gene.length > quantile(sample2gene.length, probs=.9) -> cat.b
		table(cat.a, cat.b) -> mm
		mm = mm[2:1, 2:1]
		if(mm[1,1] >= 5){
			fisher.test(mm) -> fit
			cancer.drug.tml.mat = rbind(cancer.drug.tml.mat, cbind(cancer, drug = colnames(new.ccle)[kdrug], p=fit$p.value, OR=fit$estimate) )
		} else {
			
		}
	}
}
write.table(cancer.drug.tml.mat, file="CCLE.cancer.drug.tml.FET.txt", row.names=F, quote=F, sep="\t")


