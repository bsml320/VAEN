setwd("/path/to/VAEN/Figure/Figure2")

##############################################################################
library(MASS)
library(glmnet)
library(ggplot2)
library(magrittr)

one.drugs.match = read.table("../../DATA/drugs.match.txt", as.is=T)
two.drugs.match = read.table("../../DATA/drugs.match-2.txt", as.is=T, sep="\t")

ccle.anno = read.csv("../../DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is=T)
gdsc.anno = read.delim("../../DATA/GDSC/v17.3_fitted_dose_response.txt", as.is=T)
cell.line.anno = read.csv("../../DATA/CCLE/DepMap-2018q3-celllines.csv", as.is=T)

PPs = read.table(paste("../../result/1.CCLE.latent.tsv", sep=""))
original.ss.PP = rownames(PPs)
sapply(original.ss.PP, function(x){
	new.u = u = strsplit(x, split="\\.")[[1]][1]
	if(grepl("^X", u)){
		substr(u, 2, nchar(u)) -> new.u
	}
	new.u
}) -> ss.PP
names(ss.PP) = NULL

original.ss.PP = rownames(PPs)
sapply(original.ss.PP, function(x){
	new.u = u = strsplit(x, split="\\.")[[1]]
	paste(u[3], u[4], sep="-") -> new.u
	new.u
}) -> ss.ACH
names(ss.ACH) = NULL


ccle.pred = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_CCLE.txt", header=T, as.is=T)
gdsc.pred = read.table("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_GDSC.txt", header=T, as.is=T, sep="\t")


data4plot.list = list()
shared.drugs.ori.cor = c()
for(k in 1:nrow(two.drugs.match)){
	ccle.anno.1 = ccle.anno[which(ccle.anno$Compound==two.drugs.match[k, 1]),]
	gdsc.anno.1 = gdsc.anno[which(gdsc.anno$DRUG_NAME==two.drugs.match[k, 3]),]
	print(c(two.drugs.match[k, 1], nrow(ccle.anno.1), nrow(gdsc.anno.1)))
	
	match(gdsc.anno.1$COSMIC_ID, cell.line.anno$COSMIC_ID) -> idx
	
	gdsc.anno.1.match1 = gdsc.anno.1[which( gdsc.anno.1$COSMIC_ID %in% cell.line.anno$COSMIC_ID ), ]
	gdsc.anno.1.match2 = gdsc.anno.1[which(  !gdsc.anno.1$COSMIC_ID %in% cell.line.anno$COSMIC_ID), ]
	
	match(gdsc.anno.1.match1$COSMIC_ID, cell.line.anno$COSMIC_ID) -> idx1
	match(cell.line.anno[idx1, 1], ss.ACH) -> part1.ach.idx
	gdsc.anno.1.match1 = gdsc.anno.1.match1[which(!is.na(part1.ach.idx)),]
	match(gdsc.anno.1.match1$COSMIC_ID, cell.line.anno$COSMIC_ID) -> idx1
	match(cell.line.anno[idx1, 1], ss.ACH) -> part1.ach.idx
	gdsc.Y = -gdsc.anno.1.match1[, "LN_IC50"]
	names(gdsc.Y) = gdsc.anno.1.match1[, "CELL_LINE_NAME"]
	
	
	sapply(gdsc.anno.1.match2$CELL_LINE_NAME, function(u){
		gsub("-", "", u)
	}) -> gdsc.cell.name
	union(which(gdsc.cell.name %in% cell.line.anno$Aliases), which(gdsc.anno.1.match2$CELL_LINE_NAME %in% cell.line.anno$Aliases)) -> ii
	gdsc.anno.1.match2 = gdsc.anno.1.match2[ii, ]
	gdsc.cell.name = gdsc.cell.name[ii]
	
	idx2 = c()
	for(kk in 1:nrow(gdsc.anno.1.match2)){
		if(gdsc.anno.1.match2$CELL_LINE_NAME[kk] %in% cell.line.anno$Aliases) {
			idx2 = c(idx2, match(gdsc.anno.1.match2$CELL_LINE_NAME[kk], cell.line.anno$Aliases))
		} else {
			idx2 = c(idx2, match(gdsc.cell.name[kk], cell.line.anno$Aliases))
		}
	}
	match(cell.line.anno[idx2, 1], ss.ACH) -> part2.ach.idx
	gdsc.anno.1.match2 = gdsc.anno.1.match2[which(!is.na(part2.ach.idx)),]
	
	idx2 = c()
	for(kk in 1:nrow(gdsc.anno.1.match2)){
		if(gdsc.anno.1.match2$CELL_LINE_NAME[kk] %in% cell.line.anno$Aliases) {
			idx2 = c(idx2, match(gdsc.anno.1.match2$CELL_LINE_NAME[kk], cell.line.anno$Aliases))
		} else {
			idx2 = c(idx2, match(gsub("-", "", gdsc.anno.1.match2$CELL_LINE_NAME[kk]), cell.line.anno$Aliases))
		}
	}
	
	
	match(cell.line.anno[idx2, 1], ss.ACH) -> part2.ach.idx
	new.Y = -gdsc.anno.1.match2[, "LN_IC50"]
	names(new.Y) = gdsc.anno.1.match2[, "CELL_LINE_NAME"]
	gdsc.Y = c(gdsc.Y, new.Y)
	gdsc.train.data = PPs[c(part1.ach.idx, part2.ach.idx), ]
	
	
	intersect(ccle.anno.1[,1], ss.PP) -> shared.samples
	match(shared.samples, ccle.anno.1[,1]) -> ii
	ccle.anno.2 = ccle.anno.1[ii, ]
	ccle.Y = ccle.anno.2[, "ActArea"]
	ccle.Y.IC50 = ccle.anno.2[, "IC50..uM."]
	names(ccle.Y) = ccle.anno.2[, "Primary.Cell.Line.Name"]
	match(shared.samples, ss.PP) -> ii
	ccle.train.data = PPs[ii,]
	
	shared.samples = intersect(rownames(gdsc.train.data), rownames(ccle.train.data))
	match(shared.samples, rownames(ccle.train.data)) -> ccle.ii
	match(shared.samples, rownames(gdsc.train.data)) -> gdsc.ii
	
	shared.ccle.Y = ccle.Y[ccle.ii]
	shared.ccle.Y.IC50 = ccle.Y.IC50[ccle.ii]
	shared.gdsc.Y = gdsc.Y[gdsc.ii]
	
	ccle.pred.Y = ccle.pred[match(shared.samples, ccle.pred[,1]), two.drugs.match[k, 2]]
	gdsc.pred.Y = gdsc.pred[match(shared.samples, gdsc.pred[,1]), one.drugs.match[k, 3]]
	
	cur.list = list()
	cur.list[["shared.ccle.Y"]] = shared.ccle.Y
	cur.list[["shared.gdsc.Y"]] = shared.gdsc.Y
	cur.list[["ccle.pred.Y"]] = ccle.pred.Y
	cur.list[["gdsc.pred.Y"]] = gdsc.pred.Y
	data4plot.list[[two.drugs.match[k,1] ]] = cur.list
	shared.drugs.ori.cor = rbind(shared.drugs.ori.cor, c(two.drugs.match[k, ], cor(shared.ccle.Y, shared.gdsc.Y), cor(shared.ccle.Y.IC50, shared.gdsc.Y), length(shared.ccle.Y), cor(ccle.pred.Y, gdsc.pred.Y) ) )
}
save(shared.drugs.ori.cor, data4plot.list, file="A.best.shared.drugs.ori.RData")

###########################################################################################################################
ccle = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
gdsc = read.table("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")

cancer.types = sort(unique(ccle[,2]))
shared.drugs.pred.cor = c()
for(k1 in 1:nrow(two.drugs.match)){
	df = as.data.frame(cbind(x=ccle[, two.drugs.match[k1,2] ], y=gdsc[, one.drugs.match[k1,3]]) )
	df[,1] = as.numeric(as.character(df[,1]))
	df[,2] = as.numeric(as.character(df[,2]))
	shared.drugs.pred.cor = rbind(shared.drugs.pred.cor, c(two.drugs.match[k1, ], cor(df[,1], df[,2]), nrow(df) ) )
}

pdf("2BC.shared.drugs.pdf", width=8, height=4)
par(mfrow=c(1,2), mar=c(4,4,2,1))

plot(x=shared.drugs.ori.cor[,4], y=shared.drugs.ori.cor[,7], cex=unlist(shared.drugs.ori.cor[,6])/100, xlab="Cell line original (PCC)", ylab="Cell line prediction (PCC)", xlim=c(0.1,0.9))
text(x=shared.drugs.ori.cor[,4], y=shared.drugs.ori.cor[,7], labels=shared.drugs.ori.cor[,1], pos=3)

plot(x=shared.drugs.ori.cor[,4], y=shared.drugs.pred.cor[,4], cex=unlist(shared.drugs.ori.cor[,6])/100, xlab="Cell line original (PCC)", ylab="TCGA prediction (PCC)", xlim=c(0.1,0.9))
text(x=shared.drugs.ori.cor[,4], y=shared.drugs.pred.cor[,4], labels=shared.drugs.ori.cor[,1], pos=4)

dev.off()

write.table(shared.drugs.pred.cor, file="Figure2B.txt", row.names=F, quote=F, sep="\t")
write.table(shared.drugs.ori.cor, file="Figure2C.txt", row.names=F, quote=F, sep="\t")

##########################################################################################################################

pdf("2D.FigureS3.pdf", width=10, height=10)
for(k in 1:nrow(two.drugs.match)){
	############ 1. observed
	data4plot.list[[ two.drugs.match[k,1] ]] -> cur.list
	dat.df = as.data.frame(cbind(x=cur.list$shared.ccle.Y, y=cur.list$shared.gdsc.Y))
	p1 = ggplot(data=dat.df,aes(x,y)) + 
		 stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
		 scale_fill_continuous(low="green",high="red") +
		 geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
		 guides(alpha="none") + ggtitle(paste(   one.drugs.match[k,1], "\nObserved DR, r = ", format(cor(dat.df[,1], dat.df[,2]), digits=3), ", n = ", nrow(dat.df), sep="" )) + 
		 geom_point() + labs(color="Density",fill="Density", x="Observed ActArea, CCLE", y="Observed -LN(IC50), GDSC") +
		 theme(legend.position=c(0,1),legend.justification=c(0,1), plot.title = element_text(hjust = 0.5))
	
	#################### 2. CCLE predicted
	load(paste("../../result.EN/dr.CCLE/01/1.CCLE.model.list.RData", sep="") )
	ccle.res.list = model.list[[ two.drugs.match[k,1] ]]
	Y = ccle.res.list$Ys[,1]
	pred_Y = ccle.pred[, one.drugs.match[k, 2]]
	which(Y!=-9) -> ii
	dat.df = as.data.frame(cbind(x=Y[ii], y=pred_Y[ii]))
	p2 = ggplot(data=dat.df,aes(x,y)) + 
		 stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
		 scale_fill_continuous(low="green",high="red") +
		 geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
		 guides(alpha="none") + ggtitle(paste(  two.drugs.match[k,1], "\nCCLE, r = ", format(cor(dat.df[,1], dat.df[,2]), digits=3), ", n = ", nrow(dat.df), sep="" )) + 
		 geom_point() + labs(color="Density",fill="Density", x="Observed ActArea, CCLE", y="Predicted ActArea, CCLE") +
		 theme(legend.position=c(0,1),legend.justification=c(0,1), plot.title = element_text(hjust = 0.5))
	
	#################### 3. GDSC predicted
	load(paste("../../result.EN/dr.GDSC/01/1.GDSC.model.list.RData", sep="") )
	gdsc.res.list = model.list[[ two.drugs.match[k,3] ]]
	Y = gdsc.res.list$Ys[,1]
	pred_Y = gdsc.pred[, one.drugs.match[k, 3]]
	which(Y!=-9) -> ii
	dat.df = as.data.frame(cbind(x=Y[ii], y=pred_Y[ii]))
	p3 = ggplot(data=dat.df,aes(x,y)) + 
		 stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
		 scale_fill_continuous(low="green",high="red") +
		 geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
		 guides(alpha="none") + ggtitle(paste(  two.drugs.match[k,3], "\nGDSC, r = ", format(cor(dat.df[,1], dat.df[,2]), digits=3), ", n = ", nrow(dat.df), sep="" )) + 
		 geom_point() + labs(color="Density",fill="Density", x="Observed -LN(IC50), GDSC", y="Predicted -LN(IC50), GDSC") +
		 theme(legend.position=c(0,1),legend.justification=c(0,1), plot.title = element_text(hjust = 0.5))
	
	#################### 4. predicted
	dat.df = as.data.frame(cbind(x=cur.list$ccle.pred.Y, y=cur.list$gdsc.pred.Y))
	p4 = ggplot(data=dat.df,aes(x,y)) + 
		 stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
		 scale_fill_continuous(low="green",high="red") +
		 geom_smooth(method=lm,linetype=2,colour="red",se=F) + 
		 guides(alpha="none") + ggtitle(paste(  two.drugs.match[k,3], "\nPredicted DR, r = ", format(cor(dat.df[,1], dat.df[,2]), digits=3), ", n = ", nrow(dat.df), sep="" )) + 
		 geom_point() + labs(color="Density",fill="Density", x="Predicted ActArea, CCLE", y="Predicted -LN(IC50), GDSC") +
		 theme(legend.position=c(0,1),legend.justification=c(0,1), plot.title = element_text(hjust = 0.5))
	
	multiplot(plotlist=list(p1, p2, p3, p4), layout=matrix(c(1:4), nrow=2))
	cat(two.drugs.match[k,1], ",", sep="")
}
dev.off()


###########################
library("Hmisc")

ccle = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
gdsc = read.table("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")

cancer.types = sort(unique(ccle[,2]))

shared.drugs.mat = matrix(0, nrow=nrow(one.drugs.match), ncol=length(cancer.types) )
pdf("2E.percancer.pdf", width=21, height=15)
for(k1 in 1:nrow(one.drugs.match)){
	pp.list = list()
	for(k in 1:length(cancer.types)){
		which(ccle[,2]==cancer.types[k]) -> ii;
		
		df = as.data.frame(cbind(x=ccle[ii,  one.drugs.match[k1,2] ], y=gdsc[ii, one.drugs.match[k1,3]]) )
		df[,1] = as.numeric(as.character(df[,1]))
		df[,2] = as.numeric(as.character(df[,2]))
		
		rcorr(as.matrix(df), type="pearson") -> a
		shared.drugs.mat[k1, k] = ifelse(a$r[1,2] > 0,  a$P[1,2], 1)
		
		p1 = ggplot(df, aes(x=x, y=y)) + geom_point() + ggtitle(cancer.types[k]) + xlab("CCLE") + ylab("GDSC") +
		   theme(legend.position = "none", text = element_text(size=10), axis.text.x=element_text(color="black", size = 10),axis.text.y=element_text(color="black", size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
		   geom_smooth(method="lm",color='red',data = df, aes(x=x, y=y))
	
		pp.list[[k]] = p1
	}
	
	multiplot(plotlist=pp.list, layout=matrix(c(1:35), ncol=7, byrow=T))
	cat(k1, ".", sep="")
}
dev.off()


shared.drugs.pcc.mat = matrix(0, nrow=nrow(one.drugs.match), ncol=length(cancer.types) )
for(k1 in 1:nrow(one.drugs.match)){
	for(k in 1:length(cancer.types)){
		which(ccle[,2]==cancer.types[k]) -> ii;
		
		df = as.data.frame(cbind(x=ccle[ii,  one.drugs.match[k1,2] ], y=gdsc[ii, one.drugs.match[k1,3]]) )
		df[,1] = as.numeric(as.character(df[,1]))
		df[,2] = as.numeric(as.character(df[,2]))
		
		rcorr(as.matrix(df), type="pearson") -> a
		shared.drugs.pcc.mat[k1, k] = a$r[1,2] 
	}
	cat(k1, ".", sep="")
}

###############################################################################
### 2E

rownames(shared.drugs.mat) = one.drugs.match[, 1]
colnames(shared.drugs.mat) = cancer.types
log.shared.drugs.mat = t(-log(shared.drugs.mat + 1e-16))
new = log.shared.drugs.mat[, order(apply(log.shared.drugs.mat, 2, mean))]

write.table(shared.drugs.mat, file="Figure2E.txt", sep="\t", quote=F)

library(reshape2)
log.shared.drugs.mat = t(-log(shared.drugs.mat + 1e-16))
new = log.shared.drugs.mat[, order(apply(log.shared.drugs.mat, 2, mean))]


melt(new) -> dat
pdf("2E.pdf", width=4, height=4) 
p1 = ggplot(dat, aes(x=Var2, y=value)) + 
		     geom_boxplot(outlier.shape = NA) + guides(fill=FALSE) +
		     geom_jitter(shape=21, position = position_jitter(width = 0.4), size=1) + guides(colour=FALSE) +
		     xlab("") + ylab("-log(p)") +
		     theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8, angle=90, hjust=1),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(2,2,2,2), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
print(p1)
dev.off()

