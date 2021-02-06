setwd("/path/to/VAEN/result.EN/dr.CCLE/Analysis/10.CNV")
source("../../code/unfactor.R")
library(igraph)

anno = read.table("/path/to/VAEN/DATA/UCSC/RefSeq", sep="\t", as.is=T)

###################################################################################################
dir.create("gain")
files = dir()
files = files[grep("gain.drug.limma.mat.ts.txt", files)]

###################################################################################################
for(f in 1:length(files)){
	if(file.size(files[f]) == 0)next
	cancer = strsplit(files[f], split="\\.")[[1]][1]
	dat.mat = read.table(files[f], header=F, as.is=T, sep="\t")
	
	cnv.mat.gain = read.table(paste(cancer, ".cnv_gain.txt", sep=""), header=T, as.is=T)
	
	drugs = unique(dat.mat[,3])
	cancer.gain.mat = c()
	for(kdrug in 1:length(drugs)){
		drug.sig.mat = dat.mat[dat.mat[,3]==drugs[kdrug], ]
		
		new.sensitive.mat = sensitive.mat = drug.sig.mat[which(drug.sig.mat[,6] > 1), ]
		if(nrow(sensitive.mat) > 1){
			match(sensitive.mat[,2], anno[,13]) -> ii
			cbind(sensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")

			xhetres = xhetres[order(xhetres[,11], as.numeric(as.character(xhetres[,13]))), ]  ### 11: chromocome; 13: start
			new.sensitive.mat = unfactor(xhetres)
			
		} else if(nrow(sensitive.mat) == 1){
			match(sensitive.mat[,2], anno[,13]) -> ii
			cbind(sensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")
			new.sensitive.mat = xhetres
		}
		
		
		new.insensitive.mat = insensitive.mat = drug.sig.mat[which(drug.sig.mat[,6] < 1), ]
		if(nrow(insensitive.mat) > 1){
			match(insensitive.mat[,2], anno[,13]) -> ii
			cbind(insensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")
			if(nrow(xhetres) == 0){
				new.insensitive.mat = c()
			} else {
				xhetres = xhetres[order(xhetres[,11], as.numeric(as.character(xhetres[,13]))), ]  ### 11: chromocome; 13: start
				new.insensitive.mat = unfactor(xhetres)
			}
		} else if(nrow(insensitive.mat) == 1){
			match(insensitive.mat[,2], anno[,13]) -> ii
			cbind(insensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")
			new.insensitive.mat = xhetres
		}
		
		
		new.drug.mat = rbind(new.sensitive.mat, new.insensitive.mat)
		cancer.gain.mat = rbind(cancer.gain.mat, new.drug.mat)
	}
	write.table(cancer.gain.mat, file=paste("gain/",gsub("limma.mat", "limma.mat.full", files[f]), sep=""), row.names=F, col.names=F, quote=F, sep="\t" )
}

###################################################################################################
dir.create("loss")
files = dir()
files = files[grep("loss.drug.limma.mat.ts.txt", files)]

###################################################################################################
for(f in 1:length(files)){
	if(file.size(files[f]) == 0)next
	cancer = strsplit(files[f], split="\\.")[[1]][1]
	dat.mat = read.table(files[f], header=F, as.is=T, sep="\t")
	
	cnv.mat.gain = read.table(paste(cancer, ".cnv_loss.txt", sep=""), header=T, as.is=T)
	
	drugs = unique(dat.mat[,3])
	cancer.gain.mat = c()
	for(kdrug in 1:length(drugs)){
		drug.sig.mat = dat.mat[dat.mat[,3]==drugs[kdrug], ]
		
		new.sensitive.mat = sensitive.mat = drug.sig.mat[which(drug.sig.mat[,6] > 1), ]
		if(nrow(sensitive.mat) > 1){
			match(sensitive.mat[,2], anno[,13]) -> ii
			cbind(sensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")

			xhetres = xhetres[order(xhetres[,11], as.numeric(as.character(xhetres[,13]))), ]  ### 11: chromocome; 13: start
			new.sensitive.mat = unfactor(xhetres)
			
		} else if(nrow(sensitive.mat) == 1){
			match(sensitive.mat[,2], anno[,13]) -> ii
			cbind(sensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")
			new.sensitive.mat = xhetres
		}
		
		
		new.insensitive.mat = insensitive.mat = drug.sig.mat[which(drug.sig.mat[,6] < 1), ]
		if(nrow(insensitive.mat) > 1){
			match(insensitive.mat[,2], anno[,13]) -> ii
			cbind(insensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")
			if(nrow(xhetres) == 0){
				new.insensitive.mat = c()
			} else {
				xhetres = xhetres[order(xhetres[,11], as.numeric(as.character(xhetres[,13]))), ]  ### 11: chromocome; 13: start
				new.insensitive.mat = unfactor(xhetres)
			}
		} else if(nrow(insensitive.mat) == 1){
			match(insensitive.mat[,2], anno[,13]) -> ii
			cbind(insensitive.mat, anno[ii, c(1:9, 13)]) -> xhetres
			xhetres = xhetres[!is.na(xhetres[,9]),]
			colnames(xhetres) = paste("X", 1:ncol(xhetres), sep="")
			new.insensitive.mat = xhetres
		}
		
		
		new.drug.mat = rbind(new.sensitive.mat, new.insensitive.mat)
		cancer.gain.mat = rbind(cancer.gain.mat, new.drug.mat)
	}
	write.table(cancer.gain.mat, file=paste("loss/",gsub("limma.mat", "limma.mat.full", files[f]), sep=""), row.names=F, col.names=F, quote=F, sep="\t" )
}
