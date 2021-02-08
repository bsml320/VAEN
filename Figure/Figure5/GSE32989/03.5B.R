#setwd("/path/to/GitHub/Figure/Figure5/GSE32989")
library(gplots)
library("ggplot2")

EMT = read.table("EMT.txt", as.is=T)
EMT = unique(EMT[,1])

load("GSE32989.RData")
match(EMT, rownames(expr.mat)) -> ii
ii = ii[!is.na(ii)]

prcomp(t(expr.mat[ii, ]))-> fit

which(fit$x[,1] < median(fit$x[,1])) -> idx
men.cells = rownames(fit$x)[idx]
epi.cells = rownames(fit$x)[-idx]
save(epi.cells, men.cells, file="cells.RData")

boxplot(expr.mat["ZEB1", epi.cells], expr.mat["ZEB1", men.cells])

######################

pdf("GSE32989.heatmap.pdf", width=5, height=5)
heatmap.2(expr.mat[ii, ], trace="none", col=greenred(100), cexRow=.5, cexCol=.5) -> h
dev.off()

cutree(as.hclust(h$colDendrogram),2) -> c1

########### ZEB1
dat = rbind( cbind( ZEB1=expr.mat["ZEB1", which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=expr.mat["ZEB1", which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p1 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GSE32989\n", "p = ", format(pvalue, digits=3)), x="", y = "ZEB1 expression") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1) )

#######

ccle = read.delim("CCLE.A.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p2 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("CCLE, A-model\n", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1))


#######
ccle = read.delim("CCLE.S.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p3 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("CCLE, S-model\n", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1))


############ GDSC

ccle = read.delim("GDSC.A.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p4 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GDSC, A-model\n", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1))

#######
ccle = read.delim("GDSC.S.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p5 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GDSC, S-model\n", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1))

####
source("../../../code/multiplot.R")

pdf("GSE32989.EMT.pdf", width=10, height=5)
multiplot(plotlist=list(p1,p2,p3,p4,p5), layout=matrix(1:5, nrow=1))
dev.off()


######################## Figure 7
#######
give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

ccle = read.delim("CCLE.A.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value
print(table(dat[,2]))

p2 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GSE32989, CCLE\nA-model, ", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
     stat_summary(fun.data = give.n, geom = "text")


#######
ccle = read.delim("CCLE.S.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value
print(table(dat[,2]))

p3 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GSE32989, CCLE\nS-model, ", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
     stat_summary(fun.data = give.n, geom = "text")


############ GDSC

ccle = read.delim("GDSC.A.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p4 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GSE32989, GDSC\nA-model, ", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
     stat_summary(fun.data = give.n, geom = "text")

#######
ccle = read.delim("GDSC.S.pred_GSE32989.txt", as.is=T)

dat = rbind( cbind( ZEB1=ccle$Erlotinib[which(c1==2)], grp="Epithelial-like" ),
             cbind( ZEB1=ccle$Erlotinib[which(c1==1)], grp="Mesenchymal-like" )  )

dat = as.data.frame(dat)
dat[,1] = as.numeric(as.character(dat[,1]))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p5 = ggplot(dat, aes(x=grp, y=ZEB1, fill=grp)) + geom_boxplot() + 
     labs(title=paste("GSE32989, GDSC\nS-model, ", "p = ", format(pvalue, digits=3)), x="", y = "Response to Erlotinib") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
     stat_summary(fun.data = give.n, geom = "text")

####
source("../../../code/multiplot.R")

pdf("5B.GSE32989.EMT.pdf", width=8.5, height=5)
multiplot(plotlist=list(p2,p3,p4,p5), layout=matrix(1:4, nrow=1))
dev.off()
