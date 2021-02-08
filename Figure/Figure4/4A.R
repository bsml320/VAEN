setwd("/path/to/VAEN/Figure/Figure4")

give.n <- function(x){
   return(c(y = max(x), label = length(x)))
}

###################################################################################################

drug.ccle = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss
sample.type = substr(drug.ccle[,1], 14, 15)

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	tmp.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, tmp.ccle)
}
drug.ccle = cancer.drug.ccle

###################################################################################################

########## CLIN
clin.data = read.delim("BRCA_clinicalMatrix", header=T, sep="\t")
type = substr(clin.data[,1], 14, 15)
clin.data = clin.data[which(type==type.code), ]
c(1, grep("RFS", colnames(clin.data)), grep("OS", colnames(clin.data))) -> ii
clin.data2 = clin.data[, ii]
clin.ss = gsub("\\.","-",clin.data2[,1])
rownames(clin.data2) = clin.ss
	
##############
	
which(clin.data[, "lab_proc_her2_neu_immunohistochemistry_receptor_status"] %in% c("Equivocal", "Negative", "Positive") ) -> ii
stat = clin.data[ii, c("sampleID", "lab_proc_her2_neu_immunohistochemistry_receptor_status")]
stat = stat[which(stat[,1] %in% drug.ccle[,1]), ]
	
brca.ccle = drug.ccle[match(stat[,1], drug.ccle[,1]), ]
X = stat[,2]

dat4plot.ccle = as.data.frame(cbind(Lapatinib= as.numeric( brca.ccle[,"Lapatinib"] ), X=X))
dat4plot.ccle[,1] = as.numeric(as.character(dat4plot.ccle[,1]))

dat4plot.ccle[,2] <- factor(dat4plot.ccle[,2], levels = c("Negative", "Equivocal", "Positive"))

summary(glm(dat4plot.ccle[,1] ~ as.numeric(dat4plot.ccle[,2]))) -> sfit
p = coef(sfit)[2, 4]

#> print(sfit)
#
#Call:
#glm(formula = dat4plot.ccle[, 1] ~ as.numeric(dat4plot.ccle[, 
#    2]))
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-0.98708  -0.11637   0.00865   0.12506   0.61310  
#
#Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                    1.486556   0.014168 104.927   <2e-16 ***
#as.numeric(dat4plot.ccle[, 2]) 0.074649   0.008121   9.192   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.03634106)
#
#    Null deviance: 35.886  on 904  degrees of freedom
#Residual deviance: 32.816  on 903  degrees of freedom
#AIC: -427.62
#
#Number of Fisher Scoring iterations: 2

p1 = ggplot(dat4plot.ccle, aes(x=X, y=Lapatinib, fill=X)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) + 
  labs(title=paste("CCLE-based Model\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Lapatinib") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")


###################################################################################################
###################################################################################################

drug.gdsc = read.delim("../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T)
colnames(drug.gdsc)[3:ncol(drug.gdsc)] -> drugs
cancer.types = unique(drug.gdsc[,2])
ss = gsub("\\.", "-", drug.gdsc[,1])
drug.gdsc[,1] = ss
sample.type = substr(drug.gdsc[,1], 14, 15)

cancer.drug.gdsc = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	tmp.gdsc = drug.gdsc[which(drug.gdsc[,2] == cancer & sample.type == type.code), ]
	cancer.drug.gdsc = rbind(cancer.drug.gdsc, tmp.gdsc)
}
drug.gdsc = cancer.drug.gdsc

###################################################################################################

########## CLIN
clin.data = read.delim("BRCA_clinicalMatrix", header=T, sep="\t")
type = substr(clin.data[,1], 14, 15)
clin.data = clin.data[which(type==type.code), ]
c(1, grep("RFS", colnames(clin.data)), grep("OS", colnames(clin.data))) -> ii
clin.data2 = clin.data[, ii]
clin.ss = gsub("\\.","-",clin.data2[,1])
rownames(clin.data2) = clin.ss

##############

which(clin.data[, "lab_proc_her2_neu_immunohistochemistry_receptor_status"] %in% c("Equivocal", "Negative", "Positive") ) -> ii
stat = clin.data[ii, c("sampleID", "lab_proc_her2_neu_immunohistochemistry_receptor_status")]
stat = stat[which(stat[,1] %in% drug.gdsc[,1]), ]
	
brca.gdsc = drug.gdsc[match(stat[,1], drug.gdsc[,1]), ]
X = stat[,2]

dat4plot.gdsc = as.data.frame(cbind(Lapatinib= as.numeric( brca.gdsc[,"Lapatinib"] ), X=X))
dat4plot.gdsc[,1] = as.numeric(as.character(dat4plot.gdsc[,1]))

dat4plot.gdsc[,2] <- factor(dat4plot.gdsc[,2], levels = c("Negative", "Equivocal", "Positive"))

summary(glm(dat4plot.gdsc[,1] ~ as.numeric(dat4plot.gdsc[,2]))) -> sfit
p = coef(sfit)[2, 4]


print(sfit)
#> print(sfit)
#
#Call:
#glm(formula = dat4plot.gdsc[, 1] ~ as.numeric(dat4plot.gdsc[, 
#    2]))
#
#Deviance Residuals: 
#     Min        1Q    Median        3Q       Max  
#-1.57146  -0.21490   0.04283   0.25637   1.03240  
#
#Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                    -1.68183    0.02786 -60.375  < 2e-16 ***
#as.numeric(dat4plot.gdsc[, 2])  0.09614    0.01597   6.021 2.52e-09 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for gaussian family taken to be 0.1404964)
#
#    Null deviance: 131.96  on 904  degrees of freedom
#Residual deviance: 126.87  on 903  degrees of freedom
#AIC: 796.15
#
#Number of Fisher Scoring iterations: 2

p2 = ggplot(dat4plot.gdsc, aes(x=X, y=Lapatinib, fill=X)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) +
  labs(title=paste("GDSC-based Model\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Lapatinib") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
  stat_summary(fun.data = give.n, geom = "text")


source("../../code/multiplot.R")
pdf("4A.ERBB2.Lapatinib.pdf", height=4, width=5)
multiplot(plotlist=list(p1, p2), layout=matrix(c(1,2), nrow=1))
dev.off()
