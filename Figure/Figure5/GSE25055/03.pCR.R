setwd("/path/to/VAEN/Figure/Figure5/GSE25055")
library("survminer")
library("survival")
source("../../code/multiplot.R")

load("GSE25055.RData")

######################################################################
ccle = read.table("CCLE.A.pred_GSE25055.txt", as.is=T, header=T)

t.test(ccle[ii,"Paclitaxel"] ~ pheno.anno[ii, "characteristics_ch1.18"])

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.18"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[which(dat[,2]==1),2] = "Rx Insensitive"
dat[which(dat[,2]==2),2] = "Rx Sensitive"

#dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p1 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("CCLE\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))

############################################################################
ccle = read.delim("GDSC.A.pred_GSE25055.txt", as.is=T)

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.18"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[which(dat[,2]==1),2] = "Rx Insensitive"
dat[which(dat[,2]==2),2] = "Rx Sensitive"

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p2 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GDSC\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))


pdf("CCLE.Paclitaxel.sensitive.ggplot.pdf", width=6, height=6)
multiplot(plotlist=list(p1,p2), layout=matrix(1:2, nrow=1))
dev.off()

#############################################################

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))
dat[which(dat[,2]==1),2] = "pCR"
dat[which(dat[,2]==2),2] = "RD"

dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p1 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("All Samples\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))

###
which(pheno.anno[,"characteristics_ch1.3"]=="er_status_ihc: N") -> ii
dat = data.frame(cbind(Response = ccle[ii, "Paclitaxel"], pCR=pheno.anno[ii,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[which(dat[,2]==1),2] = "pCR"
dat[which(dat[,2]==2),2] = "RD"
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p2 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("ER Negative\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))

###
which(pheno.anno[,"characteristics_ch1.3"]=="er_status_ihc: P" & pheno.anno[,"characteristics_ch1.21"] != "NA") -> ii
dat = data.frame(cbind(Response = ccle[ii, "Paclitaxel"], pCR=pheno.anno[ii,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[which(dat[,2]==1),2] = "pCR"
dat[which(dat[,2]==2),2] = "RD"
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p3 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("ER Positive\n", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))


pdf("CCLE.Paclitaxel.pCR.ggplot.pdf", width=6, height=6)
multiplot(plotlist=list(p1,p3,p2), layout=matrix(1:3, nrow=1))
dev.off()


#################################################

ccle = read.delim("A.F1-W5-PCC.GDSC.best.pred_GSE25055.txt", as.is=T)

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))
dat[which(dat[,2]==1),2] = "pCR"
dat[which(dat[,2]==2),2] = "RD"

dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p1 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("All Samples\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))

###
which(pheno.anno[,"characteristics_ch1.3"]=="er_status_ihc: N") -> ii
dat = data.frame(cbind(Response = ccle[ii, "Paclitaxel"], pCR=pheno.anno[ii,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[which(dat[,2]==1),2] = "pCR"
dat[which(dat[,2]==2),2] = "RD"
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p2 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("ER Negative\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))

###
which(pheno.anno[,"characteristics_ch1.3"]=="er_status_ihc: P") -> ii
dat = data.frame(cbind(Response = ccle[ii, "Paclitaxel"], pCR=pheno.anno[ii,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[which(dat[,2]==1),2] = "pCR"
dat[which(dat[,2]==2),2] = "RD"
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p3 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("ER Positive\n", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5))

source("C:/Users/pjia/UTH/code/multiplot.R")
multiplot(plotlist=list(p1,p3,p2), layout=matrix(1:3, nrow=1))

pdf("GDSC.Paclitaxel.pCR.ggplot.pdf", width=6, height=6)
multiplot(plotlist=list(p1,p3,p2), layout=matrix(1:3, nrow=1))
dev.off()
