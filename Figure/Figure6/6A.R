#setwd("/path/to/VAEN/Figure/Figure6")
library("gplots")

#############################################################################################################

load("../../result.EN/dr.CCLE/01/1.CCLE.model.list.RData")
ccle.model.list = model.list

obsd.ccle.mat = c()
for(k in 1:length(ccle.model.list)){
	res.list = ccle.model.list[[k]]
	drug = names(ccle.model.list)[k]
	Ys = res.list$Ys
	
	obsd.ccle.mat = cbind(obsd.ccle.mat, Ys[,1])
}
colnames(obsd.ccle.mat) = names(ccle.model.list)

apply(obsd.ccle.mat, 1, function(u)sum(u==-9)) -> check
obsd.ccle.mat = obsd.ccle.mat[check < 1, ]

### Plot 1
x = as.matrix(obsd.ccle.mat)
apply(x, 2, scale) -> x1
heatmap.2(x1, trace="none", col=bluered(75) ) -> h2.2

#############################################################################################################
### predicted CCLE

self.pred = read.table(  paste("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.txt", sep=""), as.is=T, header=T, sep="\t")
dim(self.pred)
colnames(self.pred) -> dd
dd[dd=="X17.AAG"] = "17-AAG"
colnames(self.pred) = gsub("\\.", "-", dd)

apply(self.pred[, 2:ncol(self.pred)], 1, function(u)sum(u==-9)) -> check
self.pred = self.pred[check<1, ]

x = as.matrix(self.pred[, 2:25])
apply(x, 2, scale) -> x1
heatmap.2(x1, trace="none", col=bluered(75)) -> predCCLE.hm.2

#############################################################################################################
### imputed CCLE

self.pred = read.table(  paste("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.full.txt", sep=""), as.is=T, header=T, sep="\t")
dim(self.pred)
colnames(self.pred) -> dd
dd[dd=="X17.AAG"] = "17-AAG"
colnames(self.pred) = gsub("\\.", "-", dd)

apply(self.pred[, 2:ncol(self.pred)], 1, function(u)sum(u==-9)) -> check
self.pred = self.pred[check<1, ]

x = as.matrix(self.pred[, 2:25])
apply(x, 2, scale) -> x1
heatmap.2(x1, trace="none", col=bluered(75)) -> imputed.CCLE.hm.2

#############################################################################################################
### predicted TCGA

drug.ccle = read.table(file="../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle) -> dd
dd[dd=="X17.AAG"] = "17-AAG"
colnames(drug.ccle) = gsub("\\.", "-", dd)

cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, blca.ccle)
}
ccle = cancer.drug.ccle


cc = read.table("../../DATA/TCGA.color.txt", as.is=T, sep="\t", comment.char="")
rowBar = rep("", nrow(ccle))
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	rowBar[which(ccle[,2]==cancer)] = cc[which(cc[,1]==cancer),2]
}

x = as.matrix(ccle[, 3:26])
apply(x, 2, scale) -> x1
heatmap.2(x1, trace="none", col=bluered(75), RowSideColors=rowBar) -> ccle.h.2

####################################################################################################################################

cc = rep("black", 24)
cc[which(colnames(obsd.ccle.mat)  %in% c("Lapatinib", "ZD-6474", "Erlotinib", "AZD0530") )] = "orange"
cc[which(colnames(obsd.ccle.mat)  %in% c("PD.0325901", "PD-0325901", "AZD6244") )] = "red"
cc[which(colnames(obsd.ccle.mat)  %in% c("RAF265", "PLX4720") )] = "pink"
cc[which(colnames(obsd.ccle.mat)  %in% c("PF2341066", "PHA-665752") )] = "cyan"
cc[which(colnames(obsd.ccle.mat)  %in% c("PD-0332991", "Nutlin-3") )] = "green"
cc[which(colnames(obsd.ccle.mat)  %in% c("Irinotecan", "Topotecan", "Paclitaxel") )] = "blue"
drug.cc = cc

library(ape)
pdf("6A.pdf", width=4, height=6)
par(mfrow=c(4,1), mar=c(1,5,2,5))
plot(as.phylo(as.hclust(h2.2$colDendrogram)), tip.color=drug.cc, direction="d", font=1, label.offset=0, cex=1, srt=-180, adj=1)
mtext("Cell lines, observed", cex=.8)

plot(as.phylo(as.hclust(predCCLE.hm.2$colDendrogram)), tip.color=drug.cc, direction="d", font=1, label.offset=0.5, cex=1, srt=-180, adj=1)
mtext("Cell lines, predicted", cex=.8)

plot(as.phylo(as.hclust(imputed.CCLE.hm.2$colDendrogram)), tip.color=drug.cc, direction="d", font=1, label.offset=0.5, cex=1, srt=-180, adj=1)
mtext("Cell lines, imputed", cex=.8)

plot(as.phylo(as.hclust(ccle.h.2$colDendrogram)), tip.color=drug.cc, direction="d", font=1, label.offset=0.5, cex=1, srt=-180, adj=1)
mtext("TCGA, predicted", cex=.8)

segments(1,1,3.5,1,lwd=2)
text(2.5,0.5,"EGFRi")

segments(5,1,6.5,1,lwd=2)
text(6.5,0.5,"MEKi")

segments(15,1,16.5,1,lwd=2)
text(16.5,0.5,"BRAFi")

segments(9,1,11.8,1,lwd=2)
text(10,0.5,"Cytotoxic")

segments(12,1,14,1,lwd=2)
text(13,0.5,"c-METi")

dev.off()


