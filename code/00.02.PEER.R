setwd("/work/")
library(peer)

TCGA.mat = read.table("TCGA.4VAE.tsv", header=T, as.is=T)
CCLE.mat = read.table("CCLE.4VAE.tsv", header=T, as.is=T)
VALI.mat = read.table("GSE65185.4VAE.tsv", header=T, as.is=T)

comb.mat = rbind(CCLE.mat, TCGA.mat, VALI.mat)

############ use the 2000 most variably expressed genes for PEER 
apply(comb.mat, 2, var) -> geneVar

geneVar = sort(geneVar, decreasing=T)
genes = names(geneVar)[1:2000]
expr = comb.mat[, genes]
dim(expr)

model = PEER()
PEER_setPhenoMean(model,as.matrix(expr))
dim(PEER_getPhenoMean(model))
PEER_setNk(model,10)
PEER_getNk(model)
PEER_update(model)
factors = PEER_getX(model)
dim(factor)
dim(factors)
weights = PEER_getW(model)
dim(weights)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

tiff("peer.tif")
plot(factors[,1], factors[,2])
points(factors[1:1100,1], factors[1:1100,2], col="red")
points(factors[1101:11559,1], factors[1101:11559,2], col="blue")
points(factors[11560:11629,1], factors[11560:11629,2], col="cyan")
dev.off()


######################### regress out the first 3 principal patterns
resid.comb.mat = matrix(0, nrow=nrow(comb.mat), ncol=ncol(comb.mat))
dimnames(resid.comb.mat) = dimnames(comb.mat)
for(k in 1:ncol(comb.mat)){
	summary(lm(comb.mat[,k] ~ factors[,1] + factors[,2] + factors[,3])) -> sfit
	resid.comb.mat[, k] = residuals(sfit)
	if(k %% 1000 == 0)cat(k, ".")
}

##################################################################
resid.CCLE.mat = resid.comb.mat[1:1100, ]
zeroone.resid.CCLE.mat = apply(resid.CCLE.mat, 2, minmax_normalization)
dimnames(zeroone.resid.CCLE.mat) = dimnames(CCLE.mat)
write.table(zeroone.resid.CCLE.mat, file=paste("CCLE.4VAE-peer.tsv", sep=""), row.names=T, quote=F, sep="\t")

resid.TCGA.mat = resid.comb.mat[1101:11559, ]
zeroone.resid.TCGA.mat = apply(resid.TCGA.mat, 2, minmax_normalization)
dimnames(zeroone.resid.TCGA.mat) = dimnames(TCGA.mat)
write.table(zeroone.resid.TCGA.mat, file=paste("TCGA.4VAE-peer.tsv", sep=""), row.names=T, quote=F, sep="\t")

resid.VALI.mat = resid.comb.mat[11560:11629, ]
zeroone.resid.VALI.mat = apply(resid.VALI.mat, 2, minmax_normalization)
dimnames(zeroone.resid.VALI.mat) = dimnames(VALI.mat)
write.table(zeroone.resid.VALI.mat, file=paste("GSE65185.4VAE-peer.tsv", sep=""), row.names=T, quote=F, sep="\t")

################################################################
