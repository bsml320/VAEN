setwd("/path/to/VAEN/Figure/Figure5/GSE32646")

### key genes
cur.genes = read.table("../../../result/key.genes.txt", as.is=T)
cur.genes = cur.genes[,1]

load("GSE32646.RData")
val.RPKM = expr.mat

match(cur.genes, rownames(val.RPKM)) -> ii
summary(ii)

#> summary(ii)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      5    5061   11324   10604   16771   21652     302

new.RPKM = val.RPKM[ii, ]
apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes


### scale to z, per-sample
rank.GSE32646.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE32646.gene.mat = apply(rank.GSE32646.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	u
})

scaled.GSE32646.gene.mat = apply(rank.GSE32646.gene.mat, 2, function(u){qnorm(u)} )

##########################################################################
write.table(t(scaled.GSE32646.gene.mat), file=paste("GSE32646.RANK.T.tsv", sep=""), row.names=T, quote=F, sep="\t")
write.table(scaled.GSE32646.gene.mat, file=paste("GSE32646.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################

give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

##########################################################################
##########################################################################

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

GSE32646.pred = read.table("GSE32646.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE32646.probabilities = predict(fit, as.matrix(GSE32646.pred[,-1]), s = 'lambda.min')
GSE32646.pred.mat = cbind(GSE32646.pred[,1], GSE32646.probabilities)
colnames(GSE32646.pred.mat) = c("Sample", drug)

write.table(GSE32646.pred.mat, file=paste("CCLE.A.pred_GSE32646.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

### plot
load("GSE32646.RData")
ccle = read.table("CCLE.A.pred_GSE32646.txt", header=T, as.is=T)
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,50]))
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p3 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GSE32646, CCLE\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")

print(p3)

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

GSE32646.pred = read.table("GSE32646.RANK.latent.tsv", header=T, sep="\t", as.is=T)
GSE32646.probabilities = predict(fit, as.matrix(GSE32646.pred[,-1]), s = 'lambda.min')
GSE32646.pred.mat = cbind(GSE32646.pred[,1], GSE32646.probabilities)
colnames(GSE32646.pred.mat) = c("Sample", drug)

write.table(GSE32646.pred.mat, file=paste("GDSC.A.pred_GSE32646.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

### plot
gdsc = read.table("GDSC.A.pred_GSE32646.txt", header=T, as.is=T)
dat = data.frame(cbind(Response = gdsc[, "Paclitaxel"], pCR=pheno.anno[,50]))
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p4 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GSE32646, GDSC\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")

print(p4)

pdf("5G.2.CCLE.Paclitaxel.sensitive.ggplot.pdf", width=6, height=6)
multiplot(plotlist=list(p3,p4), layout=matrix(1:2, nrow=1))
dev.off()

