setwd("/work/Figures/Figure5")
library("RColorBrewer")
library("reshape2")
library(gplots)

########################
### RUN /work/result.EN/dr.CCLE/04-mix/12.TML/04.12.CCLE.TumorMutationLoad.R
########################

dat.mat = read.table("/work/result.EN/dr.CCLE/04-mix/12.TML/CCLE.cancer.drug.tml.wilcox.txt", header=T, as.is=T)

dat.mat = dat.mat[which(dat.mat[,1]!="BRCA"),]
dat.mat = dat.mat[which(dat.mat[,1]!="THCA"),]
dat.mat = dat.mat[which(dat.mat[,1]!="LGG"),]

dat.mat[which(dat.mat[,2]=="X17.AAG"), 2] = "17-AAG"
gsub("\\.", "-", dat.mat[,2]) -> ss
dat.mat[,2] = ss


dat.mat$logP = -log10(dat.mat[,3]) * ifelse(dat.mat[,4] > 1, 1, -1)
dcast(dat.mat, cancer~drug, value.var="logP") -> mm
rownames(mm) = mm[, 1]
mm = mm[, -1]
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(mm),col=colfunc(50),trace="none")

pdf("Figure5B.CCLE.TML.mix.pdf")
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(mm),col=colfunc(20),trace="none")
rect(0,0,1,1)
dev.off()

####################################################################
