setwd("/path/to/VAEN/Figure/Figure4")

###################################################################################################
drug.ccle = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
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

library("survival")
library("survminer")

response = read.delim("../../DATA/response/drug_response.txt", as.is=T)
response = response[  which(response$drug.name == "Paclitaxel" & response$cancers == "BRCA"), ]

match(response[,2], substr(drug.ccle[,1], 1, 12)) -> ii
cbind(response, drug.ccle[ii, ] ) -> new2
new2 = new2[!is.na(ii), ]
dim(new2)


brca.clin.data = read.delim("BRCA_clinicalMatrix", as.is=T)
match(new2[,2], substr(brca.clin.data[,1], 1, 12)) -> ii
brca.clin.data = brca.clin.data[ii, ]

drug = "Paclitaxel"

samples = brca.clin.data[,1]
match(samples, drug.ccle[,1]) -> ii
drug.response = drug.ccle[ii, drug]

surv.data = brca.clin.data[match(samples, brca.clin.data[,1]),]

new3 = cbind(drug.response, surv.data)
Y1 = Surv(new3[,"X_OS"], new3[,"X_OS_IND"])

xvector = ifelse(new3[,1] > median(new3[,1]), "HR", "LR")
table(xvector)

dat = data.frame(cbind(new3,xvector) )
fit = survfit( Surv(X_OS, X_OS_IND) ~ xvector, data = dat)
coxph(Y1 ~ xvector)


pdf("4C.CCLE.Paclitaxel.BRCA.pdf", width=5, height=5)

fit = survfit( Surv(as.numeric(X_OS), as.numeric(X_OS_IND)) ~ xvector, data = dat)
g1 = ggsurvplot(fit, data=dat , risk.table = TRUE,pval = TRUE,ggtheme = theme_minimal())
print(g1)

dev.off()

###################################################################################################

drug.ccle = read.table(file="../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")

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
drug.ccle = cancer.drug.ccle


response = read.delim("../../DATA/response/drug_response.txt", as.is=T)
response = response[which(response$drug.name == "Fluorouracil"), ]

match(response[,2], substr(drug.ccle[,1], 1, 12)) -> ii
cbind(response, drug.ccle[ii, ] ) -> new2
new2 = new2[!is.na(ii), ]
dim(new2)

new2 = new2[which(new2[,1]=="STAD"),]

brca.clin.data = read.delim("STAD_clinicalMatrix", as.is=T)
match(new2[,2], substr(brca.clin.data[,1], 1, 12)) -> ii
stad.clin.data = brca.clin.data[ii, ]

library(survival)

drug = "X5.Fluorouracil"

samples = stad.clin.data[,1]
match(samples, drug.ccle[,1]) -> ii

drug.response = drug.ccle[ii, drug]
surv.data = stad.clin.data[match(samples, stad.clin.data[,1]),]

new3 = cbind(drug.response, surv.data)
Y1 = Surv(new3[,"X_OS"], new3[,"X_OS_IND"])

xvector = ifelse(new3[,1] > median(new3[,1]), "HR", "LR")
table(xvector)

dat = data.frame(cbind(new3,xvector) )
fit = survfit( Surv(X_OS, X_OS_IND) ~ xvector, data = dat)
coxph(Y1 ~ xvector)

library("survminer")

pdf("4B.pdf", width=5, height=5)
g1 = ggsurvplot(fit, data=dat , risk.table = TRUE,pval = TRUE,break.time.by = 500, ggtheme = theme_minimal())
print(g1)
dev.off()
