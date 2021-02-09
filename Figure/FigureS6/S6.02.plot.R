setwd("/path/to/VAEN/Figure/FigureS6")
library(ggplot2)
library(ggpubr)

######################################################################################
load("dr.BRAF.RData")
write.table(dr.BRAF, file="dr.BRAF.txt", row.names=F, quote=F, sep="\t")
dr.BRAF = read.table("dr.BRAF.txt", as.is=T, header=T)
dr.BRAF$mutation = ifelse(dr.BRAF$mutation == 0, "WT", "MT")

dr.BRAF$grp -> group
group[which(group == "obsd")] = "Set 1"
group[which(group == "pred")] = "Set 2"
group[which(group == "full")] = "Set 3"
dr.BRAF$grp = group

dr.BRAF[which(dr.BRAF[,1]=="X17.AAG"),1] = "17-AAG"
gsub("\\.", "-", dr.BRAF[,1]) -> ss
dr.BRAF[,1] = ss

dat = dr.BRAF[dr.BRAF[,1] %in% c("AZD6244", "PD-0325901", "PLX4720", "RAF265"), ]

give.n <- function(x){
   return(c(y = min(x) , label = length(x)))
}

pp.list = list()
for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265")){
	dat = dr.BRAF[dr.BRAF[,1] == drug, ]
	dat$mutation = as.factor(dat$mutation)
	
	if(drug == "AZD6244"){
		p = ggplot(dat, aes(x=grp, y=dr, fill=mutation)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.5, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
		    ggtitle( paste(drug, ", BRAF", sep="") ) + xlab("") + ylab("ActArea") + guides(fill=FALSE) + 
		    theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
		p1 = p + stat_compare_means(aes(group = mutation), method = "t.test", label = "p.format", label.x = 0.8, size=3) +
		     stat_summary(fun.data = give.n, geom = "text", size=3)
	} else {
		p = ggplot(dat, aes(x=grp, y=dr, fill=mutation)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.5, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
		    ggtitle( paste(drug, ", BRAF", sep="") ) + xlab("") + ylab("ActArea") + guides(fill=FALSE) + 
		    theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
		p1 = p + stat_compare_means(aes(group = mutation), method = "t.test", label = "p.format", label.x = 0.8, size=3) + 
		     stat_summary(fun.data = give.n, geom = "text", size=3)
	}
	pp.list[[paste("BRAF", drug, sep="_")]] = p1
}
rm(list=setdiff(ls(), c("pp.list", give.n)))
######################################################################################

load("dr.NRAS.RData")
write.table(dr.BRAF, file="dr.NRAS.txt", row.names=F, quote=F, sep="\t")
dr.NRAS = read.table("dr.NRAS.txt", as.is=T, header=T)
dr.NRAS$mutation = ifelse(dr.NRAS$mutation == 0, "WT", "MT")

give.n <- function(x){
   return(c(y = min(x) , label = length(x)))
}

dr.NRAS$grp -> group
group[which(group == "obsd")] = "Set 1"
group[which(group == "pred")] = "Set 2"
group[which(group == "full")] = "Set 3"
dr.NRAS$grp = group

dr.NRAS[which(dr.NRAS[,1]=="X17.AAG"),1] = "17-AAG"
gsub("\\.", "-", dr.NRAS[,1]) -> ss
dr.NRAS[,1] = ss

dat = dr.NRAS[dr.NRAS[,1] %in% c("AZD6244", "PD-0325901", "PLX4720", "RAF265"), ]

for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265")){
	dat = dr.NRAS[dr.NRAS[,1] == drug, ]
	dat$mutation = as.factor(dat$mutation)

	p = ggplot(dat, aes(x=grp, y=dr, fill=mutation)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.5, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
	    ggtitle( paste(drug, ", NRAS", sep="") ) + xlab("") + ylab("ActArea") + guides(fill=FALSE) +
	    theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
	
	p1 = p + stat_compare_means(aes(group = mutation), method = "t.test", label = "p.format", label.x = 0.8, size=3) +
	     stat_summary(fun.data = give.n, geom = "text", size=3)
	pp.list[[paste("NRAS", drug, sep="_")]] = p1
}

rm(list=setdiff(ls(), c("pp.list", give.n)))
######################################################################################

load("dr.KRAS.RData")
write.table(dr.BRAF, file="dr.KRAS.txt", row.names=F, quote=F, sep="\t")
dr.KRAS = read.table("dr.KRAS.txt", as.is=T, header=T)
dr.KRAS$mutation = ifelse(dr.KRAS$mutation == 0, "WT", "MT")

give.n <- function(x){
   return(c(y = min(x) , label = length(x)))
}

dr.KRAS$grp -> group
group[which(group == "obsd")] = "Set 1"
group[which(group == "pred")] = "Set 2"
group[which(group == "full")] = "Set 3"
dr.KRAS$grp = group

dr.KRAS[which(dr.KRAS[,1]=="X17.AAG"),1] = "17-AAG"
gsub("\\.", "-", dr.KRAS[,1]) -> ss
dr.KRAS[,1] = ss

dat = dr.KRAS[dr.KRAS[,1] %in% c("AZD6244", "PD-0325901", "PLX4720", "RAF265"), ]

for(drug in c("AZD6244", "PD-0325901", "PLX4720", "RAF265")){
	dat = dr.KRAS[dr.KRAS[,1] == drug, ]
	dat$mutation = as.factor(dat$mutation)

	p = ggplot(dat, aes(x=grp, y=dr, fill=mutation)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=1, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
	    ggtitle( paste(drug, ", KRAS", sep="") ) + xlab("") + ylab("ActArea") + guides(fill=FALSE) +
	    theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
	p1 = p + stat_compare_means(aes(group = mutation), method = "t.test", label = "p.format", label.x = 0.8, size=3) +
	     stat_summary(fun.data = give.n, geom = "text", size=3)
	pp.list[[paste("KRAS", drug, sep="_")]] = p1
}

######################################################################################
source("../../code/multiplot.R")

pdf("FigureS6.pdf", width=12, height=9)
layout <- matrix(c(1:12), nrow = 3, byrow = TRUE)
multiplot(plotlist = pp.list, layout = layout)
dev.off()

######################################################################################
