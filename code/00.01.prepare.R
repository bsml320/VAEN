setwd("/work/")

### transcriptome data

ccle = read.table("/DATA/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", skip=2, sep="\t", header=T, as.is=T)
original.TCGA.RPKM = read.delim("/DATA/TCGA/ACC/HiSeqV2", as.is=T)
val.RPKM = read.table("/DATA/GSE65185/GSE65185_CuffnormFPKM.txt", header=T, as.is=T)

genes = intersect(intersect(original.TCGA.RPKM[,1], ccle[,2]), val.RPKM[,1])
print( length(genes) )

new.ccle = ccle[match(genes, ccle[,2]),]
ccle.gene.mat = as.matrix(new.ccle[, c(-1, -2)])
rownames(ccle.gene.mat) = new.ccle[,2]

### choose genes not lowly expressed, avg.RPKM > 1 in CCLE
apply(ccle.gene.mat, 1, mean) -> rowMean
which(rowMean > 1) -> expressed
ccle.gene.mat = ccle.gene.mat[expressed, ]

####################################################
## choose tissues with >= 20 cell lines
sapply(colnames(ccle.gene.mat), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt

names(tt) = NULL
table(tt) -> t1
names(which(t1 >= 20)) -> tissues.int
which(tt %in% tissues.int) -> ii

ccle.gene.mat = ccle.gene.mat[, ii]
sapply(colnames(ccle.gene.mat), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
names(tt) = NULL

#########################################
### log2-transform CCLE data
shared.ccle.gene.mat = log2(ccle.gene.mat+1)
apply(shared.ccle.gene.mat, 1, var) -> x
shared.ccle.gene.mat = shared.ccle.gene.mat[x > 0.1, ]  ### 12082018, choose genes with variance > 0.1

original.ss.PP = colnames(shared.ccle.gene.mat)
sapply(original.ss.PP, function(x){
	new.u = u = strsplit(x, split="\\.")[[1]][1]
	if(grepl("^X", u)){
		substr(u, 2, nchar(u)) -> new.u
	}
	new.u
}) -> ss.PP
names(ss.PP) = NULL

##########################################################################
### choose genes that are mostly variably expressed

apply(shared.ccle.gene.mat, 1, var) -> rowVar
names(which(rowVar > median(rowVar))) -> mad.genes
length(mad.genes)

##########################################################################
mad5000.ccle.gene.mat = shared.ccle.gene.mat[mad.genes,]  ### gene by sample

### scale to z, per-gene
scaled.ccle.gene.mat = t(apply(mad5000.ccle.gene.mat, 1, scale))
dimnames(scaled.ccle.gene.mat) = dimnames(mad5000.ccle.gene.mat)

### zero-one
minmax_normalization = function(x){(x-min(x))/(max(x)-min(x))}

zeroone_normed.ccle.gene.mat = t(apply(scaled.ccle.gene.mat, 1, minmax_normalization))
dimnames(zeroone_normed.ccle.gene.mat) = dimnames(scaled.ccle.gene.mat)

ccle.train.mat = t(zeroone_normed.ccle.gene.mat)  ### 
print(dim(ccle.train.mat))

##################################################################### PANCAN
cancer.types = dir("/data/mshao/TCGA/")
sapply(cancer.types, nchar) -> ii
cancer.types = cancer.types[which(ii <= 4)]
cancer.types = setdiff(cancer.types, c("FPPP", "LUNG"))

	
	RPKM.mat = c()
	for(k in 1:length(cancer.types)){
		cancer = cancer.types[k]
		original.LUAD.RPKM = read.delim(paste("/data/mshao/TCGA/",cancer,"/HiSeqV2", sep=""), as.is=T)
		apply(original.LUAD.RPKM[,-1],1,sum) -> rowCheck
		non0.LUAD.RPKM = original.LUAD.RPKM[which(rowCheck!=0),]
		
		shared.LUAD.RPKM = non0.LUAD.RPKM[match(mad.genes, non0.LUAD.RPKM[,1]), -1]
		rownames(shared.LUAD.RPKM) = mad.genes
		
		### scale to z, per cancer
		scaled.LUAD.RPKM = t(apply(shared.LUAD.RPKM, 1, scale )) ## by gene
		dimnames(scaled.LUAD.RPKM) = dimnames(shared.LUAD.RPKM)
		
		### zero-one
		zeroone_normed.LUAD.RPKM = t(apply(scaled.LUAD.RPKM, 1, minmax_normalization))
		dimnames(zeroone_normed.LUAD.RPKM) = dimnames(scaled.LUAD.RPKM)
		LUAD.data = t(zeroone_normed.LUAD.RPKM)
		RPKM.mat = rbind(RPKM.mat, LUAD.data)
		cat(cancer, "\t", nrow(RPKM.mat), " ", ncol(RPKM.mat), sep="")
	}
	
	apply(RPKM.mat, 2, sum) -> x
	names(which(is.na(x))) -> genes.na
	RPKM.mat = RPKM.mat[, -which(is.na(x))]
	print(dim(RPKM.mat))
	
######################################################################################################
#### validation data set. Eventually not used in the manuscript

	apply(val.RPKM[,-1],1,sum) -> rowCheck
	non0.val.RPKM = val.RPKM[which(rowCheck!=0),]
	
	shared.LUAD.RPKM = non0.val.RPKM[match(mad.genes, non0.val.RPKM[,1]), -1]
	rownames(shared.LUAD.RPKM) = mad.genes
	shared.LUAD.RPKM = log2(shared.LUAD.RPKM + 1)
	
	### scale to z, per cancer
	scaled.LUAD.RPKM = t(apply(shared.LUAD.RPKM, 1, scale )) ## by gene
	dimnames(scaled.LUAD.RPKM) = dimnames(shared.LUAD.RPKM)
	
	### zero-one
	zeroone_normed.LUAD.RPKM = t(apply(scaled.LUAD.RPKM, 1, minmax_normalization))
	dimnames(zeroone_normed.LUAD.RPKM) = dimnames(scaled.LUAD.RPKM)
	val.mat = t(zeroone_normed.LUAD.RPKM)
	
	apply(val.mat, 2, sum) -> x
	names(which(is.na(x))) -> genes.na
	val.mat = val.mat[, -which(is.na(x))]
	print(dim(val.mat))
	
	genes2 = intersect(colnames(val.mat), colnames(RPKM.mat))
	
	write.table(val.mat[, genes2], file=paste("GSE65185.4VAE.tsv", sep=""), row.names=T, quote=F, sep="\t")
	
	RPKM.mat = RPKM.mat[, genes2]
	write.table(RPKM.mat, file=paste("TCGA.4VAE.tsv", sep=""), row.names=T, quote=F, sep="\t")
	
	ccle.train.mat = ccle.train.mat[, genes2]
	print(dim(ccle.train.mat))
	write.table(ccle.train.mat, file=paste("CCLE.4VAE.tsv", sep=""), row.names=T, quote=F, sep="\t")

