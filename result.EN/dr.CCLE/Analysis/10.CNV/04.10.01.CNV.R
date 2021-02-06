setwd("/path/to/VAEN/result.EN/dr.CCLE/Analysis/10.CNV")
source("../../code/unfactor.R")

################## prepare cancer mutations

drug.ccle = read.table(file="../../VAEN_CCLE.MIX.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss
sample.type = substr(drug.ccle[,1], 14, 15)

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
### Cluster CNV
minimum.ncount = 10

cancer.types = unique(drug.ccle[,2])
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	cancer.ccle = drug.ccle[drug.ccle[,2]==cancer, ]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	fn = paste("/path/to/VAEN/DATA/TCGA/",cancer,"/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", sep="")
	if(!file.exists(fn))next
	
	cnv.mat = read.delim(fn, header=T, as.is=T)
	rownames(cnv.mat) = as.character(cnv.mat[,1])
	cnv.mat = cnv.mat[,-1]
	type = substr(colnames(cnv.mat),14,15)
	cnv.mat = cnv.mat[,which(type==type.code)]
	cnv.ss = gsub("\\.","-",colnames(cnv.mat))
	colnames(cnv.mat) = cnv.ss
	
	fixed.ss = intersect(cancer.ccle[,1], cnv.ss)
	print(c(cancer, length(fixed.ss)))
	if(length(fixed.ss) < 10){
		print(c(cancer, " skipped, n = ", length(fixed.ss)))
		next
	}
	cnv.mat   = cnv.mat[, fixed.ss]
	
	cnv.mat.gain = t(cnv.mat)
	for(k in 1:ncol(cnv.mat.gain)){
		cnv.mat.gain[which(cnv.mat.gain[,k] <= 0 ),k] = 0;  ### CNV = 0
		cnv.mat.gain[which(cnv.mat.gain[,k] == 1 ),k] = NA; ### do not consider low amplification
		cnv.mat.gain[which(cnv.mat.gain[,k] == 2),k] = 1;  ### CNV = 2
	}
	
	cnv.mat.loss = t(cnv.mat)
	for(k in 1:ncol(cnv.mat.loss)){
		cnv.mat.loss[which(cnv.mat.loss[,k] >= 0),k] = 0;
		cnv.mat.loss[which(cnv.mat.loss[,k] == -1),k] = NA;
		cnv.mat.loss[which(cnv.mat.loss[,k] == -2),k] = 1;
	}
	
	###################################################################################################

	cnv.gain = cnv.mat.gain[fixed.ss,]
	apply(cnv.gain, 2, function(u)sum(u==1, na.rm=T) ) -> l
	if(sum(l >= 10) > 0){
		cnv.gain = cnv.gain[, which(l>=10)]
		new.cnv.gain = cbind(rownames(cnv.gain), cnv.gain)
		colnames(new.cnv.gain)[1] = "Sample"
		write.table(new.cnv.gain, file=paste(cancer,".cnv_gain.txt", sep=""), row.names=F, quote=F, sep=" ")
	}
	
	cnv.loss = cnv.mat.loss[fixed.ss,]
	apply(cnv.loss, 2, function(u)sum(u==1, na.rm=T) ) -> l
	if(sum(l >= 10) > 0){
		cnv.loss = cnv.loss[, which(l>=10)]
		new.cnv.loss = cbind(rownames(cnv.loss), cnv.loss)
		colnames(new.cnv.loss)[1] = "Sample"
		write.table(new.cnv.loss, file=paste(cancer,".cnv_loss.txt", sep=""), row.names=F, quote=F, sep=" ")
	}
	###################################################################################################
	
}

######################################################################
library(parallel)
myfun = function(kgene, candidate.genes, cnv.mat, new.ccle){
	MT.ss = cnv.mat[ which(cnv.mat[, candidate.genes[kgene] ] == 1), 1]
	match(MT.ss, blca.ccle[,1]) -> MT.ii
	
	WT.ss = cnv.mat[ which(cnv.mat[, candidate.genes[kgene]] == 0), 1]
	match(WT.ss, blca.ccle[,1]) -> WT.ii
	
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

for(ct in 1:length(cancer.types)){
	gain.drug.limma.mat = c()
	loss.drug.limma.mat = c()
	
	cancer = cancer.types[ct]
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer), ]
	apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
	
	gain.fn = paste(cancer,".cnv_gain.txt", sep="")
	if(file.exists(gain.fn)){
		cnv.gain = read.table(paste(cancer,".cnv_gain.txt", sep=""), as.is=T, header=T)
		candidate.genes = colnames(cnv.gain)[-1]
		if(length(candidate.genes) < 1)next
		
		if(length(candidate.genes) > 500){
			N = ceiling(length(candidate.genes)/30)
			print(paste(N, "steps", sep=" "))
			for(istep in 1:N){
				start = (istep-1) * 30 + 1
				end = istep * 30
				if(end > length(candidate.genes))end = length(candidate.genes)
				mclapply(start:end, myfun, candidate.genes, cnv.gain, new.ccle, mc.cores=30) -> test
				for(ks in 1:length(test)){
					kgene = (istep-1) * 30 + ks
					if(candidate.genes[kgene] == "ERBB2")print(test[[ks]]$ps.twosided)
					gain.drug.limma.mat = rbind(gain.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, test[[ks]]$ps.twosided, test[[ks]]$fcs, nMut = test[[ks]]$MT, nMutProp=test[[ks]]$MT/nrow(blca.ccle)))
				}
				cat(istep, ".", sep="")
			}
		} else {
			for(kgene in 1:length(candidate.genes) ){
				MT.ss = cnv.gain[ which(cnv.gain[, candidate.genes[kgene] ] == 1), 1]
				match(MT.ss, blca.ccle[,1]) -> MT.ii
				
				WT.ss = cnv.gain[ which(cnv.gain[, candidate.genes[kgene] ] == 0), 1]
				match(WT.ss, blca.ccle[,1]) -> WT.ii
				
				apply(new.ccle, 2, function(u){
					wilcox.test(u[MT.ii], u[WT.ii])$p.value
				}) -> ps.twosided
				
				apply(new.ccle, 2, function(u){
					mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
				}) -> fcs
				
				gain.drug.limma.mat = rbind(gain.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, ps.twosided, fcs, nMut=length(MT.ss), nMutProp=length(MT.ss)/nrow(blca.ccle) ))
			}
		}
	}
	write.table(gain.drug.limma.mat, file=paste(cancer, ".gain.drug.limma.mat.ts.txt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")
	
	loss.fn = paste(cancer,".cnv_loss.txt", sep="")
	if(file.exists(loss.fn)){
	cnv.loss = read.table(paste(cancer,".cnv_loss.txt", sep=""), as.is=T, header=T)
	apply(cnv.loss[,-1], 2, sum, na.rm=T) -> check
	cnv.loss = cnv.loss[, c(1, 1+which(check >= 10))]

	candidate.genes = colnames(cnv.loss)[-1]
	if(length(candidate.genes) < 1)next
	
	if(length(candidate.genes) > 500){
		N = ceiling(length(candidate.genes)/30)
		print(paste(N, "steps", sep=" "))
		for(istep in 1:N){
			start = (istep-1) * 30 + 1
			end = istep * 30
			if(end > length(candidate.genes))end = length(candidate.genes)
			mclapply(start:end, myfun, candidate.genes, cnv.loss, new.ccle, mc.cores=30) -> test
			for(ks in 1:length(test)){
				kgene = (istep-1) * 30 + ks
				loss.drug.limma.mat = rbind(loss.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, test[[ks]]$ps.twosided, test[[ks]]$fcs, nMut = test[[ks]]$MT, nMutProp = test[[ks]]$MT/nrow(blca.ccle)))
			}
			cat(istep, ".", sep="")
		}
	} else {
		for(kgene in 1:length(candidate.genes) ){
			MT.ss = cnv.loss[ which(cnv.loss[, candidate.genes[kgene] ] == 1), 1]
			match(MT.ss, blca.ccle[,1]) -> MT.ii
			
			WT.ss = cnv.loss[ which(cnv.loss[, candidate.genes[kgene] ] == 0), 1]
			match(WT.ss, blca.ccle[,1]) -> WT.ii
			
			apply(new.ccle, 2, function(u){
				wilcox.test(u[MT.ii], u[WT.ii])$p.value
			}) -> ps.twosided
			
			apply(new.ccle, 2, function(u){
				mean(u[MT.ii], na.rm=T)-mean(u[WT.ii], na.rm=T) 
			}) -> fcs
			
			loss.drug.limma.mat = rbind(loss.drug.limma.mat, cbind(cancer, candidate.genes[kgene], drugs, ps.twosided, fcs, nMut = length(MT.ss), nMutProp = length(MT.ss)/nrow(blca.ccle) ))
		}
	}
	}
	write.table(loss.drug.limma.mat, file=paste(cancer, ".loss.drug.limma.mat.ts.txt", sep=""), row.names=F, col.names=F, quote=F, sep="\t")
}


