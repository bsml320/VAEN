setwd("/path/to/VAEN/Figure/Figure5/GSE65185")

### transcriptome data
CCLE.latent = read.table("/data1_2/jiap/projects/18-CCLE-VAE/new/V15.2/NOPEER.RANK.Sigmoid/V15.CCLE.4VAE.RANK.tsv", header=T, as.is=T)
cur.genes = colnames(CCLE.latent)

load("GSE65185.RData")
val.RPKM = expr.mat

match(cur.genes, rownames(val.RPKM)) -> ii


#> summary(ii)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      3    4806    9911   10669   17924   23305     162


new.RPKM = val.RPKM[ii, ]
apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes


### scale to z, per-sample
rank.GSE65185.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE65185.gene.mat = apply(rank.GSE65185.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	#min(u) -> m
	#u[which(u==m)] = 1/6163
	u
})

scaled.GSE65185.gene.mat = apply(rank.GSE65185.gene.mat, 2, function(u){qnorm(u)} )
print(dim(scaled.GSE65185.gene.mat))

write.table(t(scaled.GSE65185.gene.mat), file=paste("GSE65185.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################

for(k in 1:100){
	cmd = paste("python3 GSE65185.predict_VAE.py ", k, sep="")
	system(cmd)
}


##########################################################################
##########################################################################
load("../../../result.EN/dr.CCLE/dr.CCLE.S.models.RData")
drug = "PLX4720"
res.list = dr.ccle.models[[drug]]
best.index = res.list[[ "best_index" ]]
fit <- res.list$model

#cmd = paste("python3 GSE65185.predict_VAE.py ", k, sep="")
#system(cmd)

GSE65185.pred = read.table(paste("result/", best.index, ".GSE65185.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
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
best.index

#cmd = paste("python3 GSE65185.predict_VAE.py ", k, sep="")
#system(cmd)

GSE65185.pred = read.table(paste("result/", best.index, ".GSE65185.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
GSE65185.probabilities = predict(fit, as.matrix(GSE65185.pred[,-1]), s = 'lambda.min')
GSE65185.pred.mat = cbind(GSE65185.pred[,1], GSE65185.probabilities)
colnames(GSE65185.pred.mat) = c("Sample", drug)

write.table(GSE65185.pred.mat, file=paste("GDSC.A.pred_GSE65185.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
