setwd("/path/to/VAEN/main")

### transcriptome data

ccle = read.table("../DATA/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct", skip=2, sep="\t", header=T, as.is=T)
any.TCGA.RPKM = read.delim("../DATA/TCGA/ACC/HiSeqV2", as.is=T)

genes = intersect(any.TCGA.RPKM[,1], ccle[,2])
print( length(genes) )

#> print( length(genes) )
#[1] 18217

new.ccle = ccle[match(genes, ccle[,2]),]
ccle.gene.mat = as.matrix(new.ccle[, c(-1, -2)])
rownames(ccle.gene.mat) = new.ccle[,2]

### exclude lowly-expressed genes, avg.RPKM > 1 in CCLE
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
mad5000.ccle.gene.mat = log2.ccle.gene.mat[mad.genes,]  ### gene by sample
dim(mad5000.ccle.gene.mat)

#> length(mad.genes)
#[1] 6203

#> dim(mad5000.ccle.gene.mat)
#[1] 6203 1100

##########################################################################
##########################################################################

cancer.types = dir("../DATA/TCGA/")
sapply(cancer.types, nchar) -> ii
cancer.types = cancer.types[which(ii <= 4)]
cancer.types = setdiff(cancer.types, c("FPPP", "LUNG"))


RPKM.mat = c()
cancer.type.list = list()
for(k in 1:length(cancer.types)){
	cancer = cancer.types[k]
	original.TCGA.RPKM = read.delim(paste("../DATA/TCGA/",cancer,"/HiSeqV2", sep=""), as.is=T)
	
	### exclude genes with rowSum == 0
	apply(original.TCGA.RPKM[,-1],1,sum) -> rowCheck
	non0.TCGA.RPKM = original.TCGA.RPKM[which(rowCheck!=0),]
	
	shared.TCGA.RPKM = non0.TCGA.RPKM[match(mad.genes, non0.TCGA.RPKM[,1]), -1]
	rownames(shared.TCGA.RPKM) = mad.genes
	
	apply(shared.TCGA.RPKM, 1, sum) -> check
	shared.TCGA.RPKM = shared.TCGA.RPKM[!is.na(check),]
	
	if(k==1){
		cur.genes = rownames(shared.TCGA.RPKM)
	} else {
		cur.genes = intersect(cur.genes, rownames(shared.TCGA.RPKM) )
	}
	
	t.shared.TCGA.RPKM = t(shared.TCGA.RPKM)
	RPKM.mat = rbind(RPKM.mat[, cur.genes], t.shared.TCGA.RPKM[, cur.genes] ) ### cat by samples, columns are mad.genes
	cancer.type.list[[cancer]] = colnames(shared.TCGA.RPKM)
	cat(cancer, "\t", ncol(shared.TCGA.RPKM), " ", ncol(RPKM.mat), " ", nrow(RPKM.mat), "\n", sep="")
}

#ACC     79 6203 79
#BLCA    426 6203 505
#BRCA    1218 6203 1723
#CESC    308 6203 2031
#CHOL    45 6203 2076
#COAD    329 6203 2405
#DLBC    48 6203 2453
#ESCA    196 6203 2649
#GBM     172 6203 2821
#HNSC    566 6203 3387
#KICH    91 6203 3478
#KIRC    606 6203 4084
#KIRP    323 6203 4407
#LAML    173 6203 4580
#LGG     530 6203 5110
#LIHC    423 6203 5533
#LUAD    576 6203 6109
#LUSC    553 6203 6662
#MESO    87 6203 6749
#OV      308 6203 7057
#PAAD    183 6203 7240
#PCPG    187 6203 7427
#PRAD    550 6203 7977
#READ    105 6203 8082
#SARC    265 6203 8347
#SKCM    474 6203 8821
#STAD    450 6203 9271
#TGCT    156 6203 9427
#THCA    572 6203 9999
#THYM    122 6203 10121
#UCEC    201 6203 10322
#UCS     57 6203 10379
#UVM     80 6203 10459

##########################################################################
genes2 = intersect(mad.genes, colnames(RPKM.mat))

ccle.train.mat = mad5000.ccle.gene.mat[genes2,]
print(dim(ccle.train.mat))
scaled.ccle.train.mat = apply(ccle.train.mat, 2, scale)
dimnames(scaled.ccle.train.mat) = dimnames(ccle.train.mat)

#> length(genes2)
#[1] 6163
#> print(dim(ccle.train.mat))
#[1] 1100 6163

##########################################################################

scaled.RPKM.mat = t(apply(RPKM.mat, 1, scale))
dimnames(scaled.RPKM.mat) = dimnames(RPKM.mat)

##########################################################################
### dataset for NOPEER.NO01.Sigmoid
write.table(scaled.ccle.train.mat, file=paste("V15.CCLE.4VAE.SZ.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(scaled.RPKM.mat, file=paste("V15.TCGA.4VAE.SZ.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################
