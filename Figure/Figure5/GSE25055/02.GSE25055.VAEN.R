setwd("/path/to/VAEN/result.EN/dr.CCLE/GitHub/Figure/Figure5/GSE25055")

### key genes
cur.genes = read.table("/path/to/VAEN/result/key.genes.txt", as.is=T)
cur.genes = cur.genes[,1]

load("GSE25055.RData")
val.RPKM = expr.mat

grep("^RPS", rownames(val.RPKM)) -> ii.1
grep("^RPL", rownames(val.RPKM)) -> ii.2
val.RPKM = val.RPKM[-c(ii.1, ii.2),]


match(cur.genes, rownames(val.RPKM)) -> ii
new.RPKM = val.RPKM[ii, ]
apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes


### scale to z, per-sample
rank.GSE25055.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE25055.gene.mat = apply(rank.GSE25055.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	min(u) -> m
	u[which(u==m)] = 1/6163
	u
})

scaled.GSE25055.gene.mat = apply(rank.GSE25055.gene.mat, 2, function(u){qnorm(u)} )
scaled.GSE25055.gene.mat[, which(is.na(ii))] = 0
##########################################################################
write.table(t(scaled.GSE25055.gene.mat), file=paste("GSE25055.RANK.T.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(scaled.GSE25055.gene.mat, file=paste("GSE25055.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################

##########################################################################
##########################################################################
library("Matrix")
library("glmnet")

load("../../../result.EN/dr.CCLE/dr.CCLE.A.models.RData")
drug = "Paclitaxel"
res.list = dr.ccle.models[[drug]]
fit <- res.list$model
best.index = res.list[[ "best_index" ]]
print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################
			       
GSE25055.pred = read.table("GSE25055.Rank.latent.tsv", header=T, sep="\t", as.is=T)
GSE25055.probabilities = predict(fit, as.matrix(GSE25055.pred[,-1]), s = 'lambda.min')
GSE25055.pred.mat = cbind(GSE25055.pred[,1], GSE25055.probabilities)
colnames(GSE25055.pred.mat) = c("Sample", drug)

write.table(GSE25055.pred.mat, file=paste("CCLE.A.pred_GSE25055.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################
ccle = read.table("CCLE.A.pred_GSE25055.txt", as.is=T, header=T)
load("GSE25055.RData")

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))
t.test(dat[,1] ~ dat[,2])

##########################################################################
##########################################################################

load("../../../result.EN/dr.GDSC/dr.GDSC.A.models.RData")
drug = "Paclitaxel"
res.list = dr.gdsc.models[[drug]]
fit <- res.list$model
best.index = res.list[[ "best_index" ]]
print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE25055.pred = read.table("GSE25055.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE25055.probabilities = predict(fit, as.matrix(GSE25055.pred[,-1]), s = 'lambda.min')
GSE25055.pred.mat = cbind(GSE25055.pred[,1], GSE25055.probabilities)
colnames(GSE25055.pred.mat) = c("Sample", drug)

write.table(GSE25055.pred.mat, file=paste("GDSC.A.pred_GSE25055.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
############################################################
gdsc = read.table("GDSC.A.pred_GSE25055.txt", as.is=T, header=T)

### all
dat = data.frame(cbind(Response = gdsc[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))
t.test(dat[,1] ~ dat[,2])
