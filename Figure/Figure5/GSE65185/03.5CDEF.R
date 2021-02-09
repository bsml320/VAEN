setwd("/path/to/VAEN/Figure/Figure5/GSE65185")
library("survminer")
library("survival")

########################################################################################################
########################################################################################################

ccle = read.delim("CCLE.S.pred_GSE65185.txt", as.is=T)
gsub("\\.", "-", ccle[,1]) -> ss
ccle[,1] = ss

new.GSE = ccle[grep("baseline", ccle[,1]),]

pfs.dat = read.table("PFS.txt")
paste("Pt", pfs.dat[,1],"-baseline", sep="")  -> ss
pfs.dat[,1] = ss

intersect(pfs.dat[,1], new.GSE[,1]) -> ss
pfs.dat[match(ss, pfs.dat[,1]), ] -> pfs.dat
new.GSE[match(ss, new.GSE[,1]), ] -> new.GSE
cbind(pfs.dat, new.GSE[, "PLX4720"]) -> new3

Y1 = Surv(new3[,2], rep(1, nrow(new3)) )
xvector = ifelse(new3[,3] > median(new3[,3]), "HR", "LR")
new3[,3] = as.numeric(as.character(new3[,3]))

coxph(Y1 ~ xvector)

pdf(paste("5C.CCLE.S.PLX4720.pdf", sep=""), width=5, height=5)
dat = data.frame(cbind(OS_YEAR=new3[,2], OS=rep(1, nrow(new3)), X=xvector))
dat[,1] = as.numeric(as.character(dat[,1]))
fit = survfit( Surv(as.numeric(as.character(OS_YEAR)), as.numeric(OS)) ~ X, data=dat)
p1 = ggsurvplot(fit, data=dat , risk.table = TRUE,pval = TRUE,break.time.by = 50, ggtheme = theme_minimal())
print(p1)
dev.off()

##############################################################################################################

gdsc = read.delim("GDSC.A.pred_GSE65185.txt", as.is=T)
gsub("\\.", "-", gdsc[,1]) -> ss
gdsc[,1] = ss
new.GSE = gdsc[grep("baseline", gdsc[,1]),]

pfs.dat = read.table("PFS.txt")
paste("Pt", pfs.dat[,1],"-baseline", sep="")  -> ss
pfs.dat[,1] = ss

intersect(pfs.dat[,1], new.GSE[,1]) -> ss
pfs.dat[match(ss, pfs.dat[,1]), ] -> pfs.dat
new.GSE[match(ss, new.GSE[,1]), ] -> new.GSE
cbind(pfs.dat, new.GSE[, "PLX.4720"]) -> new3

Y1 = Surv(new3[,2], rep(1, nrow(new3)) )
xvector = ifelse(new3[,3] > median(new3[,3]), "HR", "LR")
new3[,3] = as.numeric(as.character(new3[,3]))

coxph(Y1 ~ xvector)

pdf(paste("5E.GDSC.A.PLX4720.pdf", sep=""), width=5, height=5)
dat = data.frame(cbind(OS_YEAR=new3[,2], OS=rep(1, nrow(new3)), X=xvector))
dat[,1] = as.numeric(as.character(dat[,1]))
fit = survfit( Surv(as.numeric(as.character(OS_YEAR)), as.numeric(OS)) ~ X, data=dat)
g1 = ggsurvplot(fit, data=dat , risk.table = TRUE,pval = TRUE,break.time.by = 50, ggtheme = theme_minimal())
print(g1)
dev.off()
