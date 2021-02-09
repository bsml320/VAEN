#setwd("/path/to/VAEN/Figure/Figure7")
source("../../code/multiplot.R")
library("ggplot2")
library("ggrepel")

##############################################################################################################

gene.drug.limma.mat = read.table("Figure7A.txt", as.is=T, header=T)
print(table(gene.drug.limma.mat$color))

which(gene.drug.limma.mat$adjp < 0.05 & gene.drug.limma.mat$FC > 0) -> ii.1
which(gene.drug.limma.mat$adjp < 0.05 & gene.drug.limma.mat$FC < 0) -> ii.2
cc = rep(0, nrow(gene.drug.limma.mat)); cc[ii.1] = 2; cc[ii.2] = 1
gene.drug.limma.mat$color = as.factor(cc)
print(table(gene.drug.limma.mat$color))

dat.plot = data.frame(FC = gene.drug.limma.mat[,5], log10p=-log10(gene.drug.limma.mat[,4]), nMutProp = gene.drug.limma.mat[,7], color = as.factor(gene.drug.limma.mat$color) )

lab = apply(gene.drug.limma.mat, 1, function(u)paste(u[1], u[2], u[3], sep="-") )
dat.plot$lab = lab

p2 = ggplot(dat.plot, aes(x=FC, y=log10p, size=nMutProp, color=color)) + geom_point(alpha=0.2) + scale_color_manual(values=c("black", "blue", "red")) + 
     theme(legend.position = c(0.9, 0.9), plot.title = element_text(hjust = 0.5, size=8), axis.title.x = element_text(size=8) ) +
     geom_vline(xintercept = 0, linetype="dashed", color = "red", size=1) + ylab("-log10(p)") + xlab("Difference in means of MT and WT samples") +
	 guides(color = FALSE) + labs(size = "% MT")

which(dat.plot$log10p > 7.5) -> ii
new.dat.plot = dat.plot[ii, ]
p3 = p2+geom_text_repel(data=new.dat.plot, aes(label=lab), color="black", size=2)

tiff("7A.tif", width=1850, height=2050, res=300)
print(p3)
dev.off()


#> gene.drug.limma.mat[which(gene.drug.limma.mat[,2]=="ARID1A" & gene.drug.limma.mat[,3]=="Topotecan"),]
#      Cancer   Gene      Drug   p.twosided         FC nMut   nMutProp        adjp color
#335     BLCA ARID1A Topotecan 7.459900e-01 0.00401245   51 0.12592593 0.885899539     0
#15887   LUAD ARID1A Topotecan 4.506081e-01 0.13810494   13 0.02554028 1.000000000     0
#85247   STAD ARID1A Topotecan 1.117748e-05 0.37268698   30 0.07281553 0.001680347     2
#96047   UCEC ARID1A Topotecan 1.628265e-02 0.26245780   27 0.16071429 0.321308198     0
#> gene.drug.limma.mat[which(gene.drug.limma.mat[,2]=="FLT3" & gene.drug.limma.mat[,3]=="Sorafenib"),]
#      Cancer Gene      Drug   p.twosided           FC nMut   nMutProp         adjp color
#14300   LAML FLT3 Sorafenib 2.027619e-05  0.062148587   37 0.21764706 2.027619e-05     2
#58700   SKCM FLT3 Sorafenib 7.000997e-01 -0.001593346   21 0.05801105 9.942637e-01     0
