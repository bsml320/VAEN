setwd("/path/to/VAEN/Figure/FigureS7")
source("../../code/multiplot.R")

ccle.expr = read.table("../../DATA/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", skip=2, sep="\t", header=T, as.is=T)
sapply(colnames(ccle.expr)[c(-1,-2)], function(u){
	strsplit(u, split="\\.")[[1]] -> v; 
	which(v == "ACH") -> ii
	paste(v[ii:(ii+1)], collapse="-")
}) -> ccle.expr.ss
names(ccle.expr.ss) = NULL
colnames(ccle.expr)[c(-1, -2)] = ccle.expr.ss

raw.mut.mat = read.delim("../../DATA/DATA/CCLE_DepMap_18q3_maf_20180718.txt", as.is=T, sep="\t", header=T)
raw.mut.mat = raw.mut.mat[raw.mut.mat$Variant_Type!="ONP",]

silence.type = c("3'UTR", "5'Flank", "5'UTR", "De_novo_Start_OutOfFrame", "IGR", "Intron", "Nonstop_Mutation", "Silent", "Start_Codon_Del", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins")
mut.mat = raw.mut.mat[!raw.mut.mat$Variant_Classification %in% silence.type, ]
mut.mat = mut.mat[mut.mat[,2]!=0,]
genes = unique(mut.mat[,1])

ccle.pred      = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.txt", header=T, as.is=T, sep="\t")
ccle.pred.full = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.full.txt", header=T, as.is=T, sep="\t")

#######################################
load(paste("../../result.EN/dr.CCLE/01/1.CCLE.model.list.RData", sep=""))

ccle.model.list = model.list
ccle.obsd = c()
for(k in 1:length(ccle.model.list)){
	res.list = ccle.model.list[[k]]
	drug = names(ccle.model.list)[k]
	Ys = res.list$Ys
	
	ccle.obsd = cbind(ccle.obsd, Ys[,1])
}
colnames(ccle.obsd) = names(ccle.model.list)
#######################################

sapply(rownames(ccle.obsd), function(u){
	strsplit(u, split="\\.")[[1]] -> v; 
	which(v == "ACH") -> ii
	paste(v[ii:(ii+1)], collapse="-")
}) -> ccle.obsd.ss
names(ccle.obsd.ss) = NULL
rownames(ccle.obsd) = ccle.obsd.ss

obsd.mut.mat = mut.mat[mut.mat[, "Broad_ID"] %in% ccle.obsd.ss, ]

#######################################
### observed
	dr.BRAF.EGFR = c()
	gene = "BRAF"
	gene.mutation = obsd.mut.mat[obsd.mut.mat[,1]==gene, ]
	MT.ss = unique(gene.mutation[, "Broad_ID"])
	
	for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265") ){
		Y = ccle.obsd[, drug]
		names(Y) = rownames(ccle.obsd)
		Y = Y[Y!=-9]
		MT.ss = intersect(MT.ss, names(Y))
		
		EGFR.expr = t(ccle.expr[match("EGFR",ccle.expr[,2]), names(Y)])
		dr.BRAF.EGFR = rbind(dr.BRAF.EGFR, cbind(drug=drug, dr=Y, mutation=ifelse(names(Y) %in% MT.ss, 1, 0), set="obsd", EGFR=EGFR.expr ) )
	}

### predicted
	
	sapply(ccle.pred[,1], function(u){
		strsplit(u, split="\\.")[[1]] -> v; 
		which(v == "ACH") -> ii
		paste(v[ii:(ii+1)], collapse="-")
	}) -> ccle.pred.ss
	names(ccle.pred.ss) = NULL
	rownames(ccle.pred) = ccle.pred.ss
	
	colnames(ccle.pred) -> dd
	dd[which(dd=="X17.AAG")] = "17-AAG"
	colnames(ccle.pred) = gsub("\\.", "-", dd)
	
	
	pred.mut.mat = mut.mat[mut.mat[, "Broad_ID"] %in% ccle.pred.ss, ]
	gene = "BRAF"
	gene.mutation = pred.mut.mat[pred.mut.mat[,1]==gene, ]
	MT.ss = unique(gene.mutation[, "Broad_ID"])
	
	for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265") ){
		Y = ccle.pred[, drug]
		names(Y) = rownames(ccle.pred)
		Y = Y[Y!=-9]
		MT.ss = intersect(MT.ss, names(Y))
		EGFR.expr = t(ccle.expr[match("EGFR",ccle.expr[,2]), names(Y)])
		dr.BRAF.EGFR = rbind(dr.BRAF.EGFR, cbind(drug=drug, dr=Y, mutation=ifelse(names(Y) %in% MT.ss, 1, 0), set="pred", EGFR=EGFR.expr ) )
	}
	
### full

	sapply(ccle.pred.full[,1], function(u){
		strsplit(u, split="\\.")[[1]] -> v; 
		which(v == "ACH") -> ii
		paste(v[ii:(ii+1)], collapse="-")
	}) -> ccle.pred.ss
	names(ccle.pred.ss) = NULL
	rownames(ccle.pred.full) = ccle.pred.ss
	
	colnames(ccle.pred.full) -> dd
	dd[which(dd=="X17.AAG")] = "17-AAG"
	colnames(ccle.pred.full) = gsub("\\.", "-", dd)

	pred.full.mut.mat = mut.mat[mut.mat[, "Broad_ID"] %in% ccle.pred.ss, ]
	gene = "BRAF"
	gene.mutation = pred.full.mut.mat[pred.full.mut.mat[,1]==gene, ]
	MT.ss = unique(gene.mutation[, "Broad_ID"])
	
	for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265") ){
		Y = ccle.pred.full[, drug]
		names(Y) = rownames(ccle.pred.full)
		Y = Y[Y!=-9]
		MT.ss = intersect(MT.ss, names(Y))
		EGFR.expr = t(ccle.expr[match("EGFR",ccle.expr[,2]), names(Y)])
		dr.BRAF.EGFR = rbind(dr.BRAF.EGFR, cbind(drug=drug, dr=Y, mutation=ifelse(names(Y) %in% MT.ss, 1, 0), set="full", EGFR=EGFR.expr ) )
	}
	
write.table(dr.BRAF.EGFR, file="dr.BRAF.EGFR.txt", row.names=F, col.names=F, quote=F, sep="\t")

##############################################################################
give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

pp.list = list()
dr.BRAF.EGFR = read.table("dr.BRAF.EGFR.txt", as.is=T)
check = TRUE
for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265") ){
	for(kset in c("obsd", "pred", "full")){
		tmp = dr.BRAF.EGFR[dr.BRAF.EGFR[,1]==drug & dr.BRAF.EGFR[,4]==kset,]
		tmp$grp = NULL
		tmp$V3 = ifelse(tmp$V3==0, "WT", "MT")
		
		grp = rep(2, nrow(tmp))
		tmp.expr = tmp[, 5]
		grp[ which(tmp.expr > quantile(tmp.expr, probs=.75)) ] = 3
		grp[ which(tmp.expr < quantile(tmp.expr, probs=.25)) ] = 1
		tmp[, "grp"] = grp
		
		Y = as.numeric(tmp[, 2])
		X1 = as.factor(tmp[,3])
		X2 = tmp[, "grp"]
		glm(Y ~ X2 + X1 + X2:X1) -> fit2; summary(fit2) -> sfit2
		pint = coef(sfit2)[3,4]
		if(kset == "obsd")label = "Set 1"
		if(kset == "pred")label = "Set 2"
		if(kset == "full")label = "Set 3"
		
		tmp = as.data.frame(tmp)
		tmp[,2] = as.numeric(as.character(tmp[,2]))
		tmp[,3] = as.factor(tmp[,3])
		tmp[,6] = as.factor(tmp[,6])

		if(check){
			p = ggplot(tmp, aes(x=grp, y=V2, fill=V3)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
				ggtitle( paste(drug, ",", label, "\np(int) = ", format(pint, digits=3), sep="") ) + xlab("") + ylab("ActArea") + guides(fill=FALSE) + 
				stat_summary(fun.data = give.n, geom = "text", size=3) + 
				theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
			check = FALSE
		} else {
			p = ggplot(tmp, aes(x=grp, y=V2, fill=V3)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
				ggtitle( paste(drug, ",", label, "\np(int) = ", format(pint, digits=3), sep="") ) + xlab("") + ylab("ActArea") + guides(fill=FALSE) + 
				stat_summary(fun.data = give.n, geom = "text", size=3) + 
				theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
			
		}
		pp.list[[paste(drug, kset, sep="_")]] = p
	}
}

layout <- matrix(c(1:12), nrow = 3, byrow = F)

pdf("FigureS7.pdf", width=12, height=9)
multiplot(plotlist = pp.list, layout = layout)
dev.off()
