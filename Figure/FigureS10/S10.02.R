setwd("/path/to/VAEN/Figure/FigureS10")

library(reshape2)
library(ggplot2)

load("S10.data.RData")

apply(BRAF.mat, 2, sum) -> colCheck
clean.BRAF.mat = BRAF.mat[,colCheck!=max(colCheck)]
colnames(clean.BRAF.mat) = paste(colnames(clean.BRAF.mat), ", BRAF", sep="")

apply(NRAS.mat, 2, sum) -> colCheck
clean.NRAS.mat = NRAS.mat[,colCheck!=max(colCheck)]
colnames(clean.NRAS.mat) = paste(colnames(clean.NRAS.mat), ", NRAS", sep="")

apply(KRAS.mat, 2, sum) -> colCheck
clean.KRAS.mat = KRAS.mat[,colCheck!=max(colCheck)]
colnames(clean.KRAS.mat) = paste(colnames(clean.KRAS.mat), ", KRAS", sep="")

colnames(HRAS.mat) = colnames(BRAF.mat)
apply(HRAS.mat, 2, sum) -> colCheck
clean.HRAS.mat = HRAS.mat[,colCheck!=max(colCheck)]
colnames(clean.HRAS.mat) = paste(colnames(clean.HRAS.mat), ", HRAS", sep="")

cbind(clean.BRAF.mat, clean.NRAS.mat, clean.KRAS.mat, clean.HRAS.mat, EGFR.mat[, "LUAD"]) -> x
colnames(x)[ncol(x)] = "LUAD, EGFR"
x = x[, sort(colnames(x))]

dat4g <- melt(as.matrix(x) )
dat4g[,3] = as.numeric( sapply(dat4g[,3], function(u)format(as.numeric(u), digits=2)) )
cutoff = 0.05/24/33

p5 <- ggplot(data = dat4g,aes(x=Var1,y=Var2)) + 
  geom_tile(aes(fill=value),colour='white') + 
  geom_text(aes(label=value),vjust=0.5, size=2) + labs(title="CCLE, catenated") +
  scale_fill_gradientn(colours=c("red","pink","white"),
                       values=c(0,cutoff,1),
                       guide = "none", name="-log10(p)", na.value = "white") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
		axis.text.y = element_text(size = 8),
		plot.title = element_text(hjust = 0.5)
		) 

pdf("S10.pdf", width=8, height=5)
print(p5)
dev.off()
