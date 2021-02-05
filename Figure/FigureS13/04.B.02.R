#setwd("/path/to/VAEN/Figure/FigureS13")
#######################################################################################################

library(ggplot2)
library(ggrepel)
source("../../code/multiplot.R")

### plot 1
dSNP.cluster.mat = read.table("S13.B.txt", header=T, as.is=T)

dat.plot = data.frame(log2FC = (as.numeric(dSNP.cluster.mat[,6])), log10p=-log10(dSNP.cluster.mat[,5]), nMutProp = as.numeric(dSNP.cluster.mat[,9]), color=as.factor(dSNP.cluster.mat$color) )
apply(dat.plot, 1, function(u)paste(u, collapse="-")) -> check
match(unique(check), check) -> ii
dat.plot = dat.plot[ii, ]

p4 = ggplot(dat.plot, aes(x=log2FC, y=log10p, size=nMutProp, color=color)) + geom_point(alpha=0.2) + scale_color_manual(values=c("black", "blue", "red")) + 
     theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size=8) ) +
	 geom_vline(xintercept = 0, linetype="dashed", color = "red", size=1) + ggtitle("(B) Association results for mutation cluster") +
	 ylab("-log10(p)") + xlab("Difference in means of MT and WT samples")

which(dat.plot$log10p > 7.5) -> ii
new.dat.plot = dat.plot[ii, ]
labels = apply(dSNP.cluster.mat[ii, ], 1, function(u)paste(u[1], u[2], u[4], sep="-") )
new.dat.plot$lab = labels

p5 = p4+geom_text_repel(data=new.dat.plot, aes(label=lab), color="black", size=2)


tiff("FigureS13.tiff", width=2600, height=1800, res=300)
layout <- matrix(c(1,2), nrow = 1, byrow = TRUE)
multiplot(plotlist = list(p3, p5), layout = layout)
dev.off()
