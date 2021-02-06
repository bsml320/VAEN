ccle = read.table("/path/to/VAEN/Figure/Figure5/GSE20194/GDSC.A.pred_GSE20194.txt", as.is=T, header=T)
load("GSE20194.RData")

### all
dat = data.frame(cbind(Response = ccle[, "Paclitaxel"], pCR = pheno.anno[,"pcr_vs_rd:ch1"]))
dat[,1] = as.numeric(as.character(dat[,1]))

dat[,2] <- factor(dat[,2], levels = c("RD", "pCR"))

pvalue = t.test(dat[,1] ~ dat[,2])$p.value

p6 = ggplot(dat, aes(x=pCR, y=Response, fill=pCR)) + geom_boxplot() + 
     labs(title=paste("GSE20194, GDSC\n", "p = ", format(pvalue, digits=3)), x="pCR Status", y = "Predicted Response to Paclitaxel") +
	 theme(legend.position = "none", plot.title = element_text(hjust=0.5)) +
	 stat_summary(fun.data = give.n, geom = "text")

print(p6)

