setwd("/path/to/work")

### transcriptome data

ccle = read.table("/work/data/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", skip=2, sep="\t", header=T, as.is=T)
original.TCGA.RPKM = read.delim("/work/data/TCGA/ACC/HiSeqV2", as.is=T)

genes = intersect(original.TCGA.RPKM[,1], ccle[,2])
print( length(genes) )

#> print( length(genes) )
#[1] 18217

new.ccle = ccle[match(genes, ccle[,2]),]
ccle.gene.mat = as.matrix(new.ccle[, c(-1, -2)])
rownames(ccle.gene.mat) = new.ccle[,2]

### choose genes with avg.RPKM > 1 in CCLE
apply(ccle.gene.mat, 1, mean) -> rowMean
which(rowMean > 1) -> expressed
ccle.gene.mat = ccle.gene.mat[expressed, ]
dim(ccle.gene.mat)

#> dim(ccle.gene.mat)
#[1] 12406  1156

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
dim(ccle.gene.mat)

#> dim(ccle.gene.mat)
#[1] 12406  1100

#########################################
### log2-transform CCLE data
log2.ccle.gene.mat = log2(ccle.gene.mat+1)

##########################################################################
### choose genes that are most variably expressed

apply(log2.ccle.gene.mat, 1, var) -> rowVar
names(which(rowVar > median(rowVar))) -> mad.genes
length(mad.genes)
mad5000.ccle.gene.mat = log2.ccle.gene.mat[mad.genes,] 
dim(mad5000.ccle.gene.mat)

#> length(mad.genes)
#[1] 6203

#> dim(mad5000.ccle.gene.mat)
#[1] 6203 1100

##########################################################################

### shrink to 0-1, per-gene
zeroone.ccle.gene.mat = t(apply(mad5000.ccle.gene.mat, 1, function(u)u/max(u, na.rm=T)))
dimnames(zeroone.ccle.gene.mat) = dimnames(mad5000.ccle.gene.mat)

### check
summary(zeroone.ccle.gene.mat[1,])

##########################################################################
##########################################################################

cancer.types = dir("/work/data/TCGA/")
sapply(cancer.types, nchar) -> ii
cancer.types = cancer.types[which(ii <= 4)]
cancer.types = setdiff(cancer.types, c("FPPP", "LUNG"))


RPKM.mat = c()
cancer.type.list = list()
for(k in 1:length(cancer.types)){
	cancer = cancer.types[k]
	original.LUAD.RPKM = read.delim(paste("/work/data/TCGA/",cancer,"/HiSeqV2", sep=""), as.is=T)
	
	### exclude genes with rowSum == 0
	apply(original.TCGA.RPKM[,-1],1,sum) -> rowCheck
	non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
	
	shared.TCGA.RPKM = non0.TCGA.RPKM[match(mad.genes, non0.TCGA.RPKM[,1]), -1]
	rownames(shared.TCGA.RPKM) = mad.genes
	
	apply(shared.TCGA.RPKM, 1, function(u)sum(is.na(u))) -> check
	if(sum(check!=0) > 0){
		shared.TCGA.RPKM = shared.TCGA.RPKM[check==0, ]
	}
	
	### scale to z, per cancer
	scaled.TCGA.RPKM = t(apply(shared.TCGA.RPKM, 1, function(u)u/max(u, na.rm=T) )) ## by gene
	dimnames(scaled.TCGA.RPKM) = dimnames(shared.TCGA.RPKM)
	
	if(k==1){
		cur.genes = rownames(scaled.TCGA.RPKM)
	} else {
		cur.genes = intersect(cur.genes, rownames(scaled.TCGA.RPKM) )
	}
	
	t.scaled.TCGA.RPKM = t(scaled.TCGA.RPKM)
	
	RPKM.mat = rbind(RPKM.mat[, cur.genes], t.scaled.TCGA.RPKM[, cur.genes] ) ### cat by samples, columns are mad.genes
	cancer.type.list[[cancer]] = colnames(scaled.TCGA.RPKM)
	cat(cancer, "\t", ncol(scaled.TCGA.RPKM), " ", ncol(RPKM.mat), " ", nrow(RPKM.mat), "\n", sep="")
}

#> print(dim(RPKM.mat))
#[1] 10459  6163

##########################################################################
genes2 = intersect(mad.genes, colnames(RPKM.mat))

ccle.train.mat = t(zeroone.ccle.gene.mat[genes2,])
print(dim(ccle.train.mat))

#> length(genes2)
#[1] 6163
#> print(dim(ccle.train.mat))
#[1] 1100 6163

##########################################################################
### dataset for NOPEER.NO01.Sigmoid
write.table(ccle.train.mat, file=paste("/work/Z01/CCLE.4VAE.Z01.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(RPKM.mat, file=paste("/work/Z01/TCGA.4VAE.Z01.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################
