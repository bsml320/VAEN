setwd("/path/to/VAEN/Figure/Figure5/GSE20194")

### key genes
cur.genes = read.table("/path/to/VAEN/result/key.genes.txt", as.is=T)
cur.genes = cur.genes[,1]

###
load("GSE20194.RData")
val.RPKM = expr.mat

match(cur.genes, rownames(val.RPKM)) -> ii

#> summary(ii)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#      2    3081    5960    6082    9145   12546    1567


new.RPKM = val.RPKM[ii, ]
apply(new.RPKM, 2, function(u){u[is.na(u)] = 0; u} ) -> raw.RPKM
rownames(raw.RPKM) = cur.genes


rank.GSE20194.gene.mat = apply(raw.RPKM, 2, function(u)rank(u)/length(u))
rank.GSE20194.gene.mat = apply(rank.GSE20194.gene.mat, 1, function(u){
	u[which(u==1)] = 6162.5/6163
	u
})

scaled.GSE20194.gene.mat = apply(rank.GSE20194.gene.mat, 2, function(u){qnorm(u)} )
scaled.GSE20194.gene.mat[, which(is.na(ii))] = 0

write.table(scaled.GSE20194.gene.mat, file=paste("GSE20194.RANK.tsv", sep=""), row.names=T, quote=F, sep="\t")
##########################################################################


##########################################################################
##########################################################################
library("Matrix")
library("glmnet")

load("/path/to/VAEN/result.EN/dr.CCLE/dr.CCLE.A.models.RData")
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

GSE20194.pred =read.table(paste("result/", best.index, ".GSE20194.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
GSE20194.probabilities = predict(fit, as.matrix(GSE20194.pred[,-1]), s = 'lambda.min')
GSE20194.pred.mat = cbind(GSE20194.pred[,1], GSE20194.probabilities)
colnames(GSE20194.pred.mat) = c("Sample", drug)

write.table(GSE20194.pred.mat, file=paste("CCLE.A.pred_GSE20194.txt", sep=""), quote=F, sep="\t", row.names=FALSE)



ccle = read.table("CCLE.A.pred_GSE20194.txt", as.is=T, header=T)

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"pcr_vs_rd:ch1"]))
dat[,1] = as.numeric(as.character(dat[,1]))
t.test(dat[,1] ~ dat[,2])


##### plot
library("ggplot2")

give.n <- function(x){
   return(c(y = max(x), label = length(x)))
}

ccle = read.table("CCLE.A.pred_GSE20194.txt", as.is=T, header=T)
load("GSE20194.RData")

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR = pheno.anno[,"pcr_vs_rd:ch1"]))
dat[,1] = as.numeric(as.character(dat[,1]))
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p5 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GSE20194, CCLE\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")

print(p5)

##########################################################################
##########################################################################

load("/path/to/VAEN/result.EN/dr.CCLE/GitHub/result.EN/dr.GDSC/dr.GDSC.A.models.RData")
drug = "Paclitaxel"
res.list = dr.gdsc.models[[drug]]
fit <- res.list$model
best.index = res.list[[ "best_index" ]]

##########################################################################
# Go to a shell, and run VAE.prediction.py using best.index obtained above
# python3 VAE.prediction.py <best.index> <path/to/GSE20194.RANK.tsv> <D:/UTH/work/18-VAE/V15.2/NOPEER.RANK.Sigmoid/result.EN/dr.CCLE/GitHub/result/>
# After finish the above python code, a new file will be generated: GSE20194.RANK.<best.index>.latent.tsv
##########################################################################

GSE20194.pred =read.table(paste("result/", best.index, ".GSE20194.latent.tsv", sep=""), header=T, sep="\t", as.is=T)
GSE20194.probabilities = predict(fit, as.matrix(GSE20194.pred[,-1]), s = 'lambda.min')
GSE20194.pred.mat = cbind(GSE20194.pred[,1], GSE20194.probabilities)
colnames(GSE20194.pred.mat) = c("Sample", drug)

write.table(GSE20194.pred.mat, file=paste("GDSC.A.pred_GSE20194.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

##### plot
gdsc = read.table("GDSC.A.pred_GSE20194.txt", as.is=T, header=T)
load("GSE20194.RData")

### all
dat = data.frame(cbind(Response = gdsc[, "Paclitaxel"], pCR = pheno.anno[,"pcr_vs_rd:ch1"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p6 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GSE20194, GDSC\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")

print(p6)



pdf("5G.3.CCLE.Paclitaxel.sensitive.ggplot.pdf", width=6, height=6)
multiplot(plotlist=list(p5,p6), layout=matrix(1:2, nrow=1))
dev.off()
