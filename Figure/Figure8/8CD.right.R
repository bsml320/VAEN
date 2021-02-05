setwd("/path/to/VAEN/Figure/Figure8")
library("ggpubr")
library("ggplot2")

######################################################################
load("../Figure7/MC3.RData")
######################################################################
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

###################################################################### ssGSEA, SKCM

PLX4720.pp.list = list()

	cancer = "SKCM"
	mut.list = cancer.mutation.list[[ cancer ]]
	new.mut.mat.3cols = mut.list[[ "mut" ]]
	fixed.ss = mut.list[[ "ss"  ]]
	gene2WT.list = mut.list[[ "gene2WT.list" ]]
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	
	new.mut.mat.3cols = new.mut.mat.3cols[which(new.mut.mat.3cols[,2] %in% fixed.ss), ]
	which(new.mut.mat.3cols[,1] == "BRAF" & grepl("V600", new.mut.mat.3cols[,3])) -> ii
	MT.ss = unique(new.mut.mat.3cols[ii, 2])
	print( length(MT.ss) )
	
	############# EXPR
	load(paste("../../DATA/ssGSEA/",cancer,".ssgsea.result.RData", sep=""))
	gsub("\\.", "-", colnames(es.mat)) -> ss
	colnames(es.mat) = ss
	
	fixed.ss = intersect(fixed.ss, intersect(colnames(es.mat), blca.ccle[,1]))
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	es.mat = es.mat[, fixed.ss]
	sig.path = "REACTOME_DOWNREGULATION_OF_ERBB2_ERBB3_SIGNALING"

	drug = "PLX4720"
	Y = as.numeric(blca.ccle[, drug])
	names(Y) = blca.ccle[,1]
			
	MT.ss = intersect(MT.ss, names(Y))
	WT.ss = setdiff(names(Y), MT.ss)
		
	EGFR.expr = es.mat[sig.path, names(Y)]
		
	X1 = ifelse(names(Y) %in% MT.ss, "MT", "WT")
	X2 = rep(2, length(EGFR.expr))
	X2[which(EGFR.expr > quantile(EGFR.expr, probs=.75)) ] = 3
	X2[which(EGFR.expr < quantile(EGFR.expr, probs=.25)) ] = 1
	
	glm(Y ~ X2 + X1 + X2:X1) -> fit2; summary(fit2) -> sfit2
	pint = coef(sfit2)[4,4]
		
	tmp = cbind(drug = drug, cancer = cancer, Y, X1, grp=X2)
	tmp = data.frame(tmp)
	tmp[,3] = as.numeric(as.character(tmp[,3]))
	print(c(cancer, drug))
	print(table(tmp[, 4], tmp[, 5]))

	p = ggplot(tmp, aes(x=grp, y=Y, fill=X1)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
		ggtitle( paste(drug, ", ", cancer, ", p(int) = ", format(pint, digits=3), "\n", sig.path, sep="") ) + xlab("") + ylab("ActArea") + guides(fill=guide_legend()) + 
		theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
		
	p1 = p + stat_compare_means(aes(group = X1), method = "t.test", label = "p.format", label.x = 0.8, size=3) +
		 scale_x_discrete(breaks=c("1","2","3"), labels=c("Q25", "Q25_75", "Q75"))
			
	if(drug == "PLX4720" )PLX4720.pp.list[[ paste(drug, "_", cancer, sep="") ]] = p1
	
###################################################################### ssGSEA, THCA
	cancer = "THCA"

	mut.list = cancer.mutation.list[[ cancer ]]
	new.mut.mat.3cols = mut.list[[ "mut" ]]
	fixed.ss = mut.list[[ "ss"  ]]
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	
	new.mut.mat.3cols = new.mut.mat.3cols[which(new.mut.mat.3cols[,2] %in% fixed.ss), ]
	which(new.mut.mat.3cols[,1] == "BRAF" & grepl("V600", new.mut.mat.3cols[,3])) -> ii
	MT.ss = unique(new.mut.mat.3cols[ii, 2])
	print( length(MT.ss) )
	
	############# EXPR
	load(paste("../../DATA/ssGSEA/",cancer,".ssgsea.result.RData", sep=""))
	gsub("\\.", "-", colnames(es.mat)) -> ss
	colnames(es.mat) = ss
	
	fixed.ss = intersect(fixed.ss, intersect(colnames(es.mat), blca.ccle[,1]))
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	es.mat = es.mat[, fixed.ss]
	
	drug = "PLX4720"
	sig.path = "REACTOME_EGFR_DOWNREGULATION"
	print(cancer)
	Y = as.numeric(blca.ccle[, drug])
	names(Y) = blca.ccle[,1]
		
	MT.ss = intersect(MT.ss, names(Y))
	WT.ss = setdiff(names(Y), MT.ss)
		
	EGFR.expr = es.mat[sig.path, names(Y)]
		
	X1 = ifelse(names(Y) %in% MT.ss, "MT", "WT")
	X2 = rep(2, length(EGFR.expr))
	X2[which(EGFR.expr > quantile(EGFR.expr, probs=.75)) ] = 3
	X2[which(EGFR.expr < quantile(EGFR.expr, probs=.25)) ] = 1
	print(table(X2))

	glm(Y ~ X2 + X1 + X2:X1) -> fit2; summary(fit2) -> sfit2
	pint = coef(sfit2)[4,4]
		
	tmp = cbind(drug = drug, cancer = cancer, Y, X1, grp=X2)
	tmp = data.frame(tmp)
	tmp[,3] = as.numeric(as.character(tmp[,3]))
	print(c(cancer, drug))
	print(table(tmp[, 4], tmp[, 5]))

	p = ggplot(tmp, aes(x=grp, y=Y, fill=X1)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
		ggtitle( paste(drug, ", ", cancer, ", p(int) = ", format(pint, digits=3), "\n", sig.path, sep="") ) + xlab("") + ylab("ActArea") + guides(fill=guide_legend()) + 
		theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
		
	p1 = p + stat_compare_means(aes(group = X1), method = "t.test", label = "p.format", label.x = 0.8, size=3) +
	     scale_x_discrete(breaks=c("1","2","3"), labels=c("Q25", "Q25_75", "Q75"))
	
	PLX4720.pp.list[[ paste(drug, "_", cancer, "_Braf", sep="") ]] = p1
	
######################################################################

layout <- matrix(c(1:2), nrow = 1, byrow = T)
pdf("8CD.right.pdf", width=8, height=3)
multiplot(plotlist = PLX4720.pp.list, layout = layout)
dev.off()

