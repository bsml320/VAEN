setwd("/path/to/GitHub/Figure/Figure5/GSE32989")

cur.genes = read.table("../../../result/key.genes.txt", as.is=T)
cur.genes = cur.genes[,1]

load("GSE32989.RData")
val.RPKM = expr.mat

match(cur.genes, rownames(val.RPKM)) -> ii
new.RPKM = val.RPKM[ii, ]

apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes

### scale to z, per-sample
rank.GSE32989.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE32989.gene.mat = apply(rank.GSE32989.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	u
})

scaled.GSE32989.gene.mat = apply(rank.GSE32989.gene.mat, 2, function(u){qnorm(u)} )
write.table(t(scaled.GSE32989.gene.mat), file=paste("GSE32989.RANK.T.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(scaled.GSE32989.gene.mat, file=paste("GSE32989.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")			       

##########################################################################
##########################################################################

load("../../../result.EN/dr.CCLE/dr.CCLE.A.models.RData")
drug = "Erlotinib"
res.list = dr.ccle.models[[drug]]
fit <- res.list$model
best.index = res.list[[ "best_index" ]]
print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE32989.pred = read.table("GSE32989.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE32989.probabilities = predict(fit, as.matrix(GSE32989.pred[,-1]), s = 'lambda.min')
GSE32989.pred.mat = cbind(GSE32989.pred[,1], GSE32989.probabilities)
colnames(GSE32989.pred.mat) = c("Sample", drug)

write.table(GSE32989.pred.mat, file=paste("CCLE.A.pred_GSE32989.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

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

GSE32989.pred = read.table("GSE32989.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE32989.probabilities = predict(fit, as.matrix(GSE32989.pred[,-1]), s = 'lambda.min')
GSE32989.pred.mat = cbind(Sample=GSE32989.pred[,1], GSE32989.probabilities)
colnames(GSE32989.pred.mat) = c("Sample", drug)

write.table(GSE32989.pred.mat, file=paste("CCLE.S.pred_GSE32989.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

##########################################################################
##########################################################################

load("../../../result.EN/dr.GDSC/dr.GDSC.A.models.RData")
drug = "Erlotinib"
res.list = dr.gdsc.models[[drug]]
fit <- res.list$model
best.index = res.list[[ "best_index" ]]
print(best.index)
##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE32989.pred = read.table("GSE32989.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE32989.probabilities = predict(fit, as.matrix(GSE32989.pred[,-1]), s = 'lambda.min')
GSE32989.pred.mat = cbind(GSE32989.pred[,1], GSE32989.probabilities)
colnames(GSE32989.pred.mat) = c("Sample", drug)

write.table(GSE32989.pred.mat, file=paste("GDSC.A.pred_GSE32989.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

##########################################################################
##########################################################################

load("../../../result.EN/dr.GDSC/dr.GDSC.S.models.RData")
drug = "Erlotinib"
res.list = dr.gdsc.models[[drug]]
best.index = res.list[[ "best_index" ]]
fit <- res.list$model

print(best.index)
##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE32989.pred = read.table("GSE32989.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE32989.probabilities = predict(fit, as.matrix(GSE32989.pred[,-1]), s = 'lambda.min')
GSE32989.pred.mat = cbind(Sample=GSE32989.pred[,1], GSE32989.probabilities)
colnames(GSE32989.pred.mat) = c("Sample", drug)

write.table(GSE32989.pred.mat, file=paste("GDSC.S.pred_GSE32989.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
