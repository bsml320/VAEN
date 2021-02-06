setwd("/path/to/VAEN/Figure/Figure4")

###################################################################################################

drug.ccle = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, blca.ccle)
}
drug.ccle = cancer.drug.ccle

###################################################################################################

cancer.ccle = drug.ccle[which(drug.ccle[,2]=="LUAD"),]

load("LUAD.MET_amp.ss.RData")
MET.exon14.ii = match(gsub("\\.", "-", ss), cancer.ccle[,1]) 
match(LUAD.MET_amp.ss, cancer.ccle[,1]) -> MET.amp.ii
MET.WT.ii = setdiff(1:nrow(cancer.ccle), MET.amp.ii)

###
dat4plot = data.frame(cbind( Y=cancer.ccle[, "PF2341066"], grp=ifelse(cancer.ccle[,1] %in% LUAD.MET_amp.ss, "MET Gain", "Other") ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

p = t.test(dat4plot[,1] ~ dat4plot[,2])$p.value

p1 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("PF2341066, MET CNV\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to PF2341066") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")

###
dat4plot = data.frame(cbind( Y=cancer.ccle[, "PHA.665752"], grp=ifelse(cancer.ccle[,1] %in% LUAD.MET_amp.ss, "MET Gain", "Other") ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))
p = t.test(dat4plot[,1] ~ dat4plot[,2])$p.value

p2 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("PHA665752, MET CNV\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to PHA665752") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")

###
load("LUAD.MET.expr.RData")

X1 = rep("1", length(gene.expr))
X1[which(gene.expr < quantile(gene.expr, probs=.25))] = "0"
X1[which(gene.expr > quantile(gene.expr, probs=.75))] = "2"

dat4plot = data.frame(cbind( Y=cancer.ccle[, "PF2341066"], grp=X1 ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

summary(glm(dat4plot[,1] ~ as.numeric(dat4plot[,2]))) -> sfit
p = coef(sfit)[2, 4]

p3 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("PF2341066, MET expression\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to PF2341066") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) + 
  scale_x_discrete(labels = c("0" = "Q25", "1" =  "Q25_75", "2" = "Q75" ) ) +
  stat_summary(fun.data = give.n, geom = "text")


###

X2 = rep("1", length(gene.expr))
X2[which(gene.expr < quantile(gene.expr, probs=.25))] = "0"
X2[which(gene.expr > quantile(gene.expr, probs=.75))] = "2"

dat4plot = data.frame(cbind( Y=cancer.ccle[, "PHA.665752"], grp=X2 ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

summary(glm(dat4plot[,1] ~ as.numeric(dat4plot[,2]))) -> sfit
p = coef(sfit)[2, 4]

p4 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("PHA665752, MET expression\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to PHA665752") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) + 
  scale_x_discrete(labels = c("0" = "Q25", "1" =  "Q25_75", "2" = "Q75" ) ) +
  stat_summary(fun.data = give.n, geom = "text")


pdf("4E.MET.CCLE.pdf", height=6, width=6)
multiplot(plotlist=list(p1,p3,p2,p4), layout=matrix(1:4, nrow=2))
dev.off()

###################################################################################################################
drug.ccle = read.table("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, blca.ccle)
}
drug.ccle = cancer.drug.ccle

###################################################################################################################

cancer.ccle = drug.ccle[which(drug.ccle[,2]=="LUAD"),]

load("LUAD.MET_amp.ss.RData")
MET.exon14.ii = match(gsub("\\.", "-", ss), cancer.ccle[,1]) 
match(LUAD.MET_amp.ss, cancer.ccle[,1]) -> MET.amp.ii
MET.WT.ii = setdiff(1:nrow(cancer.ccle), MET.amp.ii)

###
dat4plot = data.frame(cbind( Y=cancer.ccle[, "Crizotinib"], grp=ifelse(cancer.ccle[,1] %in% LUAD.MET_amp.ss, "MET Gain", "Other") ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

p = t.test(dat4plot[,1] ~ dat4plot[,2])$p.value

p1 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("Crizotinib, MET CNV\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Crizotinib") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")

###
dat4plot = data.frame(cbind( Y=cancer.ccle[, "PHA.665752"], grp=ifelse(cancer.ccle[,1] %in% LUAD.MET_amp.ss, "MET Gain", "Other") ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))
p = t.test(dat4plot[,1] ~ dat4plot[,2])$p.value

p2 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("PHA665752, MET CNV\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to PHA665752") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")

###
dat4plot = data.frame(cbind( Y=cancer.ccle[, "Foretinib"], grp=ifelse(cancer.ccle[,1] %in% LUAD.MET_amp.ss, "MET Gain", "Other") ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))
p = t.test(dat4plot[,1] ~ dat4plot[,2])$p.value

p3 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("Foretinib, MET CNV\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Foretinib") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")

###
load("LUAD.MET.expr.RData")

X2 = rep("1", length(gene.expr))
X2[which(gene.expr < quantile(gene.expr, probs=.25))] = "0"
X2[which(gene.expr > quantile(gene.expr, probs=.75))] = "2"

dat4plot = data.frame(cbind( Y=cancer.ccle[, "Crizotinib"], grp=X2 ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

summary(glm(dat4plot[,1] ~ as.numeric(dat4plot[,2]))) -> sfit
p = coef(sfit)[2, 4]

p4 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("Crizotinib, MET expression\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Crizotinib") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) + 
  scale_x_discrete(labels = c("0" = "Q25", "1" =  "Q25_75", "2" = "Q75" ) ) +
  stat_summary(fun.data = give.n, geom = "text")


###

X2 = rep("1", length(gene.expr))
X2[which(gene.expr < quantile(gene.expr, probs=.25))] = "0"
X2[which(gene.expr > quantile(gene.expr, probs=.75))] = "2"

dat4plot = data.frame(cbind( Y=cancer.ccle[, "PHA.665752"], grp=X2 ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

summary(glm(dat4plot[,1] ~ as.numeric(dat4plot[,2]))) -> sfit
p = coef(sfit)[2, 4]

p5 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("PHA665752, MET expression\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to PHA665752") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) + 
  scale_x_discrete(labels = c("0" = "Q25", "1" =  "Q25_75", "2" = "Q75" ) ) +
  stat_summary(fun.data = give.n, geom = "text")


###
X2 = rep("1", length(gene.expr))
X2[which(gene.expr < quantile(gene.expr, probs=.25))] = "0"
X2[which(gene.expr > quantile(gene.expr, probs=.75))] = "2"

dat4plot = data.frame(cbind( Y=cancer.ccle[, "Foretinib"], grp=X2 ))
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))

summary(glm(dat4plot[,1] ~ as.numeric(dat4plot[,2]))) -> sfit
p = coef(sfit)[2, 4]

p6 = ggplot(dat4plot, aes(y=Y, x=grp)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("Foretinib, MET expression\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Foretinib") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5))+ 
  scale_x_discrete(labels = c("0" = "Q25", "1" =  "Q25_75", "2" = "Q75" ) ) +
  stat_summary(fun.data = give.n, geom = "text")


pdf("4F.MET.GDSC.pdf", height=6, width=9)
multiplot(plotlist=list(p1,p4,p2,p5,p3,p6), layout=matrix(1:6, nrow=2))
dev.off()

###################################################################################################################
