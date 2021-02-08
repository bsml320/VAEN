setwd("/path/to/GitHub/Figure/Figure5/GSE33072")

cur.genes = read.table("../../../result/key.genes.txt", as.is=T)
cur.genes = cur.genes[,1]


load("GSE33072.RData")
val.RPKM = expr.mat
match(cur.genes, rownames(val.RPKM)) -> ii

#> summary(ii)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      3    4806    9911   10669   17924   23305     162

new.RPKM = val.RPKM[ii, ]
apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes

### rank-based p, per-sample
rank.GSE33072.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE33072.gene.mat = apply(rank.GSE33072.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	u
})

### p to z
scaled.GSE33072.gene.mat = apply(rank.GSE33072.gene.mat, 2, function(u){qnorm(u)} )
write.table(t(scaled.GSE33072.gene.mat), file=paste("GSE33072.RANK.T.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(scaled.GSE33072.gene.mat, file=paste("GSE33072.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")

##########################################################################
##########################################################################

load("../../../result.EN/dr.CCLE/dr.CCLE.A.models.RData")
drug = "Erlotinib"
res.list = dr.ccle.models[[drug]]
best.index = res.list[[ "best_index" ]]
fit <- res.list$model
print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE33072.pred = read.table("GSE33072.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE33072.probabilities = predict(fit, as.matrix(GSE33072.pred[,-1]), s = 'lambda.min')
GSE33072.pred.mat = cbind(GSE33072.pred[,1], GSE33072.probabilities)
colnames(GSE33072.pred.mat) = c("Sample", drug)
write.table(GSE33072.pred.mat, file=paste("CCLE.A.pred_GSE33072.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

##########################################################################
##########################################################################
load("../../../result.EN/dr.CCLE/dr.CCLE.S.models.RData")
drug = "Erlotinib"
res.list = dr.ccle.models[[drug]]
best.index = res.list[[ "best_index" ]]
fit <- res.list$model
print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE33072.pred = read.table("GSE33072.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE33072.probabilities = predict(fit, as.matrix(GSE33072.pred[,-1]), s = 'lambda.min')
GSE33072.pred.mat = cbind(Sample=GSE33072.pred[,1], GSE33072.probabilities)
colnames(GSE33072.pred.mat) = c("Sample", drug)

write.table(GSE33072.pred.mat, file=paste("CCLE.S.pred_GSE33072.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
