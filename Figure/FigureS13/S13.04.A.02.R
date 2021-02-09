#setwd("/path/to/VAEN/Figure/FigureS13")
source("../../code/multiplot.R")
library(ggplot2)
library(ggrepel)

abc = read.table("FigureS13.A.txt", as.is=T, header=T, sep="\t")
dat.plot = data.frame(FC = as.numeric(abc[,5]), log10p=-log10(abc[,4]), nMutProp = as.numeric(abc[,7]), color = as.factor(abc$color) )

lab = apply(abc, 1, function(u)paste(u[1], u[2], u[3], sep="-") )
dat.plot$lab = lab

p2 = ggplot(dat.plot, aes(x=FC, y=log10p, size=nMutProp, color=color)) + geom_point(alpha=0.2) + scale_color_manual(values=c("black", "blue", "red")) + 
      theme(legend.position = c(0.9, 0.9), plot.title = element_text(hjust = 0.5, size=8) ) +
     geom_vline(xintercept = 0, linetype="dashed", color = "red", size=1) + ylab("-log10(p)") + xlab("Difference in means of MT and WT samples") +
	 ggtitle("(A) Association results for mutations") +
	 guides(color = FALSE) + labs(size = "% MT")
	 
which(dat.plot$log10p > 8) -> ii
new.dat.plot = dat.plot[ii, ]
p3 = p2+geom_text_repel(data=new.dat.plot, aes(label=lab), color="black", size=2)

tiff("S13.A.tif", width=1800, height=2000, res=300)
print(p3)
dev.off()
