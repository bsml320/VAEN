#setwd("/path/to/VAEN/Figure/Figure7")
source("../../code/multiplot.R")
library("ggplot2")
library("ggrepel")

##############################################################################################################

gene.drug.limma.mat = read.table("Figure7A.txt", as.is=T, header=T)
print(table(gene.drug.limma.mat$color))

which(gene.drug.limma.mat$adjp < 0.05 & gene.drug.limma.mat$FC > 0) -> ii.1
pos.mat = gene.drug.limma.mat[ii.1,]
which(gene.drug.limma.mat$adjp < 0.05 & gene.drug.limma.mat$FC < 0) -> ii.2
neg.mat = gene.drug.limma.mat[ii.2,]
cc = rep(0, nrow(gene.drug.limma.mat)); cc[ii.1] = 2; cc[ii.2] = 1
gene.drug.limma.mat$color = as.factor(cc)
print(table(gene.drug.limma.mat$color))

dat.plot = data.frame(log2FC = gene.drug.limma.mat[,5], log10p=-log10(gene.drug.limma.mat[,4]), nMutProp = gene.drug.limma.mat[,7], color = as.factor(gene.drug.limma.mat$color) )

lab = apply(gene.drug.limma.mat, 1, function(u)paste(u[1], u[2], u[3], sep="-") )
dat.plot$lab = lab

p2 = ggplot(dat.plot, aes(x=log2FC, y=log10p, size=nMutProp, color=color)) + geom_point(alpha=0.2) + scale_color_manual(values=c("black", "blue", "red")) + 
     theme(legend.position = c(0.9, 0.9), plot.title = element_text(hjust = 0.5, size=8), axis.title.x = element_text(size=8) ) +
     geom_vline(xintercept = 0, linetype="dashed", color = "red", size=1) + ylab("-log10(p)") + xlab("Difference in means of MT and WT samples") +
	 guides(color = FALSE) + labs(size = "% MT")

which(dat.plot$log10p > 7.5) -> ii
new.dat.plot = dat.plot[ii, ]
p3 = p2+geom_text_repel(data=new.dat.plot, aes(label=lab), color="black", size=2)

tiff("7A.tif", width=1850, height=2050, res=300)
print(p3)
dev.off()
