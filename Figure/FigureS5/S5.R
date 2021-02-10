setwd("/path/to/VAEN/Figure/FigureS5")

#####################################################################################

load("../../result.EN/dr.CCLE/CCLE.A.info.RData")
load("../../result.EN/dr.CCLE/CCLE.S.info.RData")

############################################################################################
drugs = colnames(all.avg_CV_R2.mat)
solid.drugs = c("Erlotinib", "AZD0530", "PLX4720", "TKI258", "ZD.6474")

pdf("FigureS5.pdf", width=8, height=10.5)
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

#####################################################################################
