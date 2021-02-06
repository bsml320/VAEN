#setwd("/path/to/VAEN/Figure/Figure6")
library("RColorBrewer")
library("reshape2")
library("gplots")

########################
dat.mat = read.table("6B.CCLE.tml.ttest.txt", header=T, as.is=T)

dat.mat[which(dat.mat[,2]=="X17.AAG"), 2] = "17-AAG"
gsub("\\.", "-", dat.mat[,2]) -> ss
dat.mat[,2] = ss


dat.mat$logP = -log10(dat.mat[,3]) * ifelse(dat.mat[,4] > 1, 1, -1)
dcast(dat.mat, cancer~drug, value.var="logP") -> mm
rownames(mm) = mm[, 1]
mm = mm[, -1]
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(mm),col=colfunc(50),trace="none")


pdf("6B.pdf")
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(mm),col=colfunc(20),trace="none")
dev.off()

####################################################################
