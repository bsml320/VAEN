setwd("/path/to/VAEN/result.EN/dr.CCLE")

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")

#####################################################################################
load("../../DATA/TCGA.ss.mat.RData")
#####################################################################################
anno = read.csv("../../DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
drugs = sort(unique(anno$Compound))
#####################################################################################

load("CCLE.A.info.RData")
load("CCLE.S.info.RData")

############################################################################################
solid.drugs = c("Erlotinib", "AZD0530", "PLX4720", "TKI258", "ZD.6474")

pdf("CCLE.MIX-F1-W5-PCC.ROC.pdf", width=8, height=10.5)
par(mfrow=c(4,3), mar=c(4,4,2,1))
for(k in 1:24){
	drug = drugs[k]
	L = min( c(all.avg_CV_R2.mat[,k], solid.avg_CV_R2.mat[,k]), na.rm=T)
	H = max( c(all.avg_CV_R2.mat[,k], solid.avg_CV_R2.mat[,k]), na.rm=T)
	Lx = min( c(all.in_sample_R2.mat[,k], solid.in_sample_R2.mat[,k]), na.rm=T)
	Hx = max( c(all.in_sample_R2.mat[,k], solid.in_sample_R2.mat[,k]), na.rm=T)
	
	### model: All
	plot(x=all.in_sample_R2.mat[,k], y=all.avg_CV_R2.mat[,k], main=, xlab="", ylab="", ylim=c(L, H), xlim=c(Lx, Hx))
	points(x=solid.in_sample_R2.mat[,k], y=solid.avg_CV_R2.mat[,k], pch=3)
	ordiellipse(rbind(cbind(all.in_sample_R2.mat[,k],all.avg_CV_R2.mat[,k]), cbind(solid.in_sample_R2.mat[,k],solid.avg_CV_R2.mat[,k])), groups=c(rep(1,100),rep(2,100)),col=c(1:2), display = "sites", kind = "sd", label = F, conf = 0.95, lty=5,lwd=0.5)
	
	which.max(all.avg_CV_R2.mat[,k]) -> idx
	points(all.in_sample_R2.mat[idx,k], all.avg_CV_R2.mat[idx,k], pch=19, col="cyan")
	
	which.max(solid.avg_CV_R2.mat[,k]) -> idx
	points(solid.in_sample_R2.mat[idx,k], solid.avg_CV_R2.mat[idx,k], pch=3, col="cyan")
	
	mtext(text="In-sample PCC", side=1, line=2.3, cex=.9)
	mtext(text="CV-R2", side=2, line=2.3, cex=.9)
	mtext(text=drugs[k], side=3, line=.6)
}
dev.off()

#> solid.drugs
#[1] "AZD0530"   "Lapatinib" "LBW242"    "PLX4720"

#####################################################################################

all.TCGA.pred.mat   = read.table("VAEN_CCLE.A.pred_TCGA.txt", as.is=T, header=T)
solid.TCGA.pred.mat = read.table("VAEN_CCLE.S.pred_TCGA.txt", as.is=T, header=T)

TCGA.pred.mat = all.TCGA.pred.mat
immune.cancer = c("LAML", "DLBC", "THYM")
which(!(TCGA.pred.mat[,2] %in% immune.cancer)) -> ii

for(k in 1:length(solid.drugs)){
	drug = solid.drugs[k]
	cat("Updated ", drug, "\n", sep="")
	TCGA.pred.mat[ii, drug] = solid.TCGA.pred.mat[ii, drug]
}
write.table(TCGA.pred.mat,  file=paste("VAEN_CCLE.MIX.pred_TCGA.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################################################
all.CCLE.pred.mat   = read.table("VAEN_CCLE.A.pred_CCLE.txt", as.is=T, header=T)
solid.CCLE.pred.mat = read.table("VAEN_CCLE.S.pred_CCLE.txt", as.is=T, header=T)

pdf("MIX-F1-W5-PCC.obsd.vs.pred.pdf", width=8, height=4)
par(mfrow=c(1,2), cex=1, mar=c(4,4,3,1))
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	
	######### all
	load( paste("01/1.CCLE.model.list.RData", sep="") )
	model.list[[ drug ]] -> res.list
	Ys = res.list$Ys
	which(Ys[,1]!=-9) -> ii
	r = cor(Ys[ii, 1], all.CCLE.pred.mat[ii, kdrug+1])
	plot(Ys[ii, 1], all.CCLE.pred.mat[ii, kdrug+1], main=paste(drug, ", all, PCC = ", format(r, digits=3), sep=""), xlab="Observed CCLE", ylab="Predicted CCLE (all)")
	
	######### solid
	load( paste("01S/1.CCLE.model.list.S.RData", sep="") )
	model.list[[ drug ]] -> res.list
	Ys = res.list$Ys
	which(Ys[,1]!=-9) -> ii
	r = cor(Ys[ii, 1], solid.CCLE.pred.mat[ii, kdrug+1])
	plot(Ys[ii, 1], solid.CCLE.pred.mat[ii, kdrug+1], main=paste(drug, ", solid, PCC = ", format(r, digits=3), sep=""), xlab="Observed CCLE", ylab="Predicted CCLE (solid)")
}
dev.off()


############################################################################################
############################################################################################

all.CCLE.pred.mat   = read.table("VAEN_CCLE.A.pred_CCLE.txt", as.is=T, header=T)
solid.CCLE.pred.mat = read.table("VAEN_CCLE.S.pred_CCLE.txt", as.is=T, header=T)

CCLE.pred.mat = all.CCLE.pred.mat
for(k in 1:length(solid.drugs)){
	drug = solid.drugs[k]
	CCLE.pred.mat[, drug] = solid.CCLE.pred.mat[, drug]
}
write.table(CCLE.pred.mat,  file=paste("VAEN_CCLE.MIX.pred_CCLE.txt", sep=""), quote=F, sep="\t", row.names=FALSE)

############################################################################################
############################################################################################


all.CCLE.pred.mat   = read.table("VAEN_CCLE.A.pred_CCLE.full.txt", as.is=T, header=T)
solid.CCLE.pred.mat = read.table("VAEN_CCLE.S.pred_CCLE.full.txt", as.is=T, header=T)

CCLE.pred.mat = all.CCLE.pred.mat
for(k in 1:length(solid.drugs)){
	drug = solid.drugs[k]
	CCLE.pred.mat[, drug] = solid.CCLE.pred.mat[, drug]
}
write.table(CCLE.pred.mat,  file=paste("VAEN_CCLE.MIX.pred_CCLE.full.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
