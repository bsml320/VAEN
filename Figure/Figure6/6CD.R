#setwd("/path/to/VAEN/Figure/Figure6")
library("RColorBrewer")
library(reshape2)
library(gplots)
library(ggplot2)

give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

##########################################################################################
drug.ccle = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")
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
##########################################################################################

pdf("6CD.pdf", width=5.5, height=4)
for(kdrug in 1:length(drugs)){
	dat4plot = c()
	for(ct in 1:length(cancer.types)){
		cancer = cancer.types[ct]
		cancer = cancer.types[ct]
		blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer), ]
		
		
		if(cancer == "LAML"){
			fn = paste("../../MC3/LAML_wustl", sep="")
		} else {
			fn = paste("../../MC3/",cancer,"_mc3.txt", sep="")
		}
		mut.mat = read.delim(fn, as.is=T)
		
		fixed.ss = intersect(blca.ccle[,1], mut.mat[,1] )
		if(length(fixed.ss) < 50){
			next
		}
		
		tapply(mut.mat[,1], mut.mat[,1], length) -> sample2gene.length
		sample2gene.length = sample2gene.length[fixed.ss]
		sample2gene.length[which(is.na(sample2gene.length))] = 0
		names(sample2gene.length) = fixed.ss
	
		blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
		apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
		
		sample2gene.length.log = log(sample2gene.length)
		
		which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.75)) -> ii
		dat4plot = rbind(dat4plot, cbind(cancer, drug=colnames(new.ccle)[kdrug], group = 0, value=sample2gene.length.log[ii] ) )
		dat4plot = rbind(dat4plot, cbind(cancer, drug=colnames(new.ccle)[kdrug], group = 1, value=sample2gene.length.log[-ii]  ) )
	}
	dat4plot = dat4plot[!is.na(dat4plot[,4]), ]
	dat4plot = as.data.frame(dat4plot)
	dat4plot[,4] = as.numeric(as.character(dat4plot[,4]))
	
	p = ggplot(dat4plot, aes(x=cancer, y=value, color=group)) + geom_boxplot( aes(color = group)) + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjus=1), plot.title = element_text(hjust = 0.5) ) + ggtitle(drugs[kdrug]) + theme(legend.position="none") + 
		ylab("log(TMB)") + xlab("") +
		stat_summary(fun.data = give.n, geom = "text", size=2)
	print(p)	
}

p = ggplot(dat4plot, aes(x=cancer, y=value, fill=group)) + geom_boxplot( aes(color = group)) 
print(p)	
dev.off()

