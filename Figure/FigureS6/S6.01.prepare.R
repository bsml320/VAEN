setwd("/path/to/VAEN/Figure/FigureS6")

raw.mut.mat = read.delim("../../DATA/CCLE/CCLE_DepMap_18q3_maf_20180718.txt", as.is=T, sep="\t", header=T)
raw.mut.mat = raw.mut.mat[raw.mut.mat$Variant_Type!="ONP",]

silence.type = c("3'UTR", "5'Flank", "5'UTR", "De_novo_Start_OutOfFrame", "IGR", "Intron", "Nonstop_Mutation", "Silent", "Start_Codon_Del", "Start_Codon_Ins", "Start_Codon_SNP", "Stop_Codon_Del", "Stop_Codon_Ins")
mut.mat = raw.mut.mat[!raw.mut.mat$Variant_Classification %in% silence.type, ]
mut.mat = mut.mat[mut.mat[,2]!=0,]
genes = unique(mut.mat[,1])

#######################################
load("../../result.EN/dr.CCLE/01/1.CCLE.model.list.RData")
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

ccle.pred      = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.txt", header=T, as.is=T, sep="\t")
ccle.pred.full = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.full.txt", header=T, as.is=T, sep="\t")

###################################################################################################
sapply(rownames(ccle.obsd), function(u){
	strsplit(u, split="\\.")[[1]] -> v; 
	which(v == "ACH") -> ii
	paste(v[ii:(ii+1)], collapse="-")
}) -> ccle.obsd.ss
names(ccle.obsd.ss) = NULL
rownames(ccle.obsd) = ccle.obsd.ss

obsd.mut.mat = mut.mat[mut.mat[, "Broad_ID"] %in% ccle.obsd.ss, ]
### observed

genes = c("EGFR", "BRAF", "KRAS", "NRAS")
for(gene in genes){
	dr.BRAF = c()
	#gene = "BRAF"
	gene.mutation = obsd.mut.mat[obsd.mut.mat[,1]==gene, ]
	MT.ss = unique(gene.mutation[, "Broad_ID"])
	
	for(kdrug in 1:ncol(ccle.obsd)){
		drug = colnames(ccle.obsd)[kdrug]
		Y = ccle.obsd[, kdrug]
		names(Y) = rownames(ccle.obsd)
		Y = Y[Y!=-9]
		MT.ss = intersect(MT.ss, names(Y))
		dr.BRAF = rbind(dr.BRAF, cbind(drug=drug, dr=Y, mutation=ifelse(names(Y) %in% MT.ss, 1, 0), grp="obsd" ) )
	}

	###################################################################################################
	sapply(ccle.pred[,1], function(u){
		strsplit(u, split="\\.")[[1]] -> v; 
		which(v == "ACH") -> ii
		paste(v[ii:(ii+1)], collapse="-")
	}) -> ccle.pred.ss
	names(ccle.pred.ss) = NULL
	rownames(ccle.pred) = ccle.pred.ss

	pred.mut.mat = mut.mat[mut.mat[, "Broad_ID"] %in% ccle.pred.ss, ]
	### observed
	gene.mutation = pred.mut.mat[pred.mut.mat[,1]==gene, ]
	MT.ss = unique(gene.mutation[, "Broad_ID"])
	if(length(MT.ss) < 5)next
	
	for(kdrug in 2:ncol(ccle.pred)){
		drug = colnames(ccle.pred)[kdrug]
		Y = ccle.pred[, kdrug]
		names(Y) = rownames(ccle.pred)
		Y = Y[Y!=-9]
		MT.ss = intersect(MT.ss, names(Y))
		dr.BRAF = rbind(dr.BRAF, cbind(drug=drug, dr=Y, mutation=ifelse(names(Y) %in% MT.ss, 1, 0), grp="pred" ) )
	}
	
	###################################################################################################
	sapply(ccle.pred.full[,1], function(u){
		strsplit(u, split="\\.")[[1]] -> v; 
		which(v == "ACH") -> ii
		paste(v[ii:(ii+1)], collapse="-")
	}) -> ccle.pred.ss
	names(ccle.pred.ss) = NULL
	rownames(ccle.pred.full) = ccle.pred.ss

	pred.full.mut.mat = mut.mat[mut.mat[, "Broad_ID"] %in% ccle.pred.ss, ]

	### observed
	gene.mutation = pred.full.mut.mat[pred.full.mut.mat[,1]==gene, ]
	MT.ss = unique(gene.mutation[, "Broad_ID"])
	
	for(kdrug in 2:ncol(ccle.pred.full)){
		drug = colnames(ccle.pred.full)[kdrug]
		Y = ccle.pred.full[, kdrug]
		names(Y) = rownames(ccle.pred.full)
		Y = Y[Y!=-9]
		MT.ss = intersect(MT.ss, names(Y))
		dr.BRAF = rbind(dr.BRAF, cbind(drug=drug, dr=Y, mutation=ifelse(names(Y) %in% MT.ss, 1, 0), grp="full" ) )
	}

	save(dr.BRAF, file=paste("dr.",gene,".RData", sep=""))
}
