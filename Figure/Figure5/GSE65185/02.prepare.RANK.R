setwd("/path/to/VAEN/Figure/Figure5/GSE65185")

### transcriptome data
cur.genes = read.table("../../../result/key.genes.txt", as.is=T)
cur.genes = cur.genes[,1]

val.RPKM = read.table("GSE65185_CuffnormFPKM.txt", header=T, as.is=T)

match(cur.genes, val.RPKM[,1]) -> ii

#> summary(ii)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      4    5276   12183   12188   19338   25266      52


new.RPKM = val.RPKM[ii, -1]
apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes

log2.raw.RPKM = log2(raw.RPKM + 1)

### scale to z, per-sample
rank.GSE65185.gene.mat = apply(log2.raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE65185.gene.mat = apply(rank.GSE65185.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	u
})

scaled.GSE65185.gene.mat = apply(rank.GSE65185.gene.mat, 2, function(u){qnorm(u)} )
print(dim(scaled.GSE65185.gene.mat))

##########################################################################
### dataset for NOPEER.NO01.Sigmoid
write.table(t(scaled.GSE65185.gene.mat), file=paste("GSE65185.RANK.T.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(scaled.GSE65185.gene.mat, file=paste("GSE65185.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################

##########################################################################

## for(k in 1:100){
##	cmd = paste("python3 GSE65185.predict_VAE.py ", k, sep="")
##	system(cmd)
## }


##########################################################################
##########################################################################
load("../../../result.EN/dr.CCLE/dr.CCLE.S.models.RData")
drug = "PLX4720"
res.list = dr.ccle.models[[drug]]
best.index = res.list[[ "best_index" ]]
fit <- res.list$model

print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE65185.pred = read.table("GSE65185.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE65185.probabilities = predict(fit, as.matrix(GSE65185.pred[,-1]), s = 'lambda.min')
GSE65185.pred.mat = cbind(Sample=GSE65185.pred[,1], GSE65185.probabilities)
colnames(GSE65185.pred.mat) = c("Sample", drug)

write.table(GSE65185.pred.mat, file=paste("CCLE.S.pred_GSE65185.txt", sep=""), quote=F, sep="\t", row.names=FALSE)



##########################################################################
##########################################################################

load("../../../result.EN/dr.GDSC/dr.GDSC.A.models.RData")
drug = "PLX-4720"
res.list = dr.gdsc.models[[drug]]
fit <- res.list$model
best.index = res.list[[ "best_index" ]]
print(best.index)

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> </path/to/VAEN/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE65185.pred = read.table("GSE65185.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE65185.probabilities = predict(fit, as.matrix(GSE65185.pred[,-1]), s = 'lambda.min')
GSE65185.pred.mat = cbind(GSE65185.pred[,1], GSE65185.probabilities)
colnames(GSE65185.pred.mat) = c("Sample", drug)

write.table(GSE65185.pred.mat, file=paste("GDSC.A.pred_GSE65185.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
