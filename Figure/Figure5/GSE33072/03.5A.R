setwd("/path/to/GitHub/Figure/Figure5/GSE33072")
library("survminer")
library("survival")

load("GSE33072.RData")
grep("progression-free survival time",pheno.anno[,"characteristics_ch1.7"]) -> ii
pheno.anno[ii, ] -> pheno.anno.pfs

########################################################################################################
ccle = read.delim("CCLE.A.pred_GSE33072.txt", as.is=T, header=T)

OS_YEAR = as.numeric(sapply(as.character(pheno.anno.pfs[,"characteristics_ch1.7"]), function(u)trimws(strsplit(u, split=":")[[1]][2])))
OS = sapply( as.character(pheno.anno.pfs[,"characteristics_ch1.6"]), function(u)trimws(strsplit(u, split=":")[[1]][2]))

match(pheno.anno.pfs[,2], ccle[,1]) -> ii
ccle = ccle[ii, ]

Y1 = Surv(as.numeric(OS_YEAR), as.numeric(OS))
Y_response = ccle[, "Erlotinib"]
xvector = ifelse(Y_response > median(Y_response), "HR", "LR")

fit = survfit( Surv(as.numeric(OS_YEAR), as.numeric(OS)) ~ xvector)
coxph(Y1 ~ xvector)


pdf(paste("5A.CCLE.A.Erlotinib.pdf", sep=""), width=5, height=5)
dat = data.frame(cbind(OS_YEAR=OS_YEAR, OS=OS, X=xvector))
dat[,1] = as.numeric(as.character(dat[,1]))
fit = survfit( Surv(as.numeric(OS_YEAR), as.numeric(OS)) ~ X, data=dat)
g1 = ggsurvplot(fit, data=dat , risk.table = TRUE,pval = TRUE,break.time.by = 1, ggtheme = theme_minimal())
print(g1)
dev.off()

########################################################################################################
ccle = read.delim("CCLE.S.pred_GSE33072.txt", as.is=T, header=T)

OS_YEAR = as.numeric(sapply(as.character(pheno.anno.pfs[,"characteristics_ch1.7"]), function(u)trimws(strsplit(u, split=":")[[1]][2])))
OS = sapply( as.character(pheno.anno.pfs[,"characteristics_ch1.6"]), function(u)trimws(strsplit(u, split=":")[[1]][2]))

match(pheno.anno.pfs[,2], ccle[,1]) -> ii
ccle = ccle[ii, ]

Y1 = Surv(as.numeric(OS_YEAR), as.numeric(OS))
Y_response = ccle[, "Erlotinib"]
xvector = ifelse(Y_response > median(Y_response), "HR", "LR")

fit = survfit( Surv(as.numeric(OS_YEAR), as.numeric(OS)) ~ xvector)
coxph(Y1 ~ xvector)


pdf(paste("5A.CCLE.S.Erlotinib.pdf", sep=""), width=5, height=5)
dat = data.frame(cbind(OS_YEAR=OS_YEAR, OS=OS, X=xvector))
dat[,1] = as.numeric(as.character(dat[,1]))
fit = survfit( Surv(as.numeric(OS_YEAR), as.numeric(OS)) ~ X, data=dat)
g2 = ggsurvplot(fit, data=dat , risk.table = TRUE,pval = TRUE,break.time.by = 1, ggtheme = theme_minimal())
print(g2)
dev.off()
