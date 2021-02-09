setwd("/path/to/VAEN/Figure/Figure5/GSE25055")
library("survminer")
library("survival")
source("../../../code/multiplot.R")
give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

load("GSE25055.RData")

######################################################################
ccle = read.table("CCLE.A.pred_GSE25055.txt", as.is=T, header=T)

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))
t.test(dat[,1] ~ dat[,2])

dat[which(dat[,2]=="dlda30_prediction: pCR"),2] = "pCR"
dat[which(dat[,2]=="dlda30_prediction: RD"),2] = "RD"
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p1 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("CCLE\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")

############################################################################
gdsc = read.delim("GDSC.A.pred_GSE25055.txt", as.is=T)

### all
dat = data.frame(cbind(Response = gdsc[, "Paclitaxel"], pCR=pheno.anno[,"characteristics_ch1.21"]))
dat[,1] = as.numeric(as.character(dat[,1]))
t.test(dat[,1] ~ dat[,2])

dat[which(dat[,2]=="dlda30_prediction: pCR"),2] = "pCR"
dat[which(dat[,2]=="dlda30_prediction: RD"),2] = "RD"
dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p2 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GDSC\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")


pdf("5G.1.CCLE.Paclitaxel.sensitive.ggplot.pdf", width=6, height=6)
multiplot(plotlist=list(p1,p2), layout=matrix(1:2, nrow=1))
dev.off()
