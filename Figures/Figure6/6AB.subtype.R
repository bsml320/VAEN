setwd("/work/Figures/Figure6/")
source("/work/code/multiplot.R")
library(ggplot2)

load("/work/data/UCSC/ncbiRefSeq.01252019.hg19.gene2ll.RData")
long.genes = names(gene2ll[gene2ll > 200000])


load("/work/result.EN/dr.CCLE/04-mix/01.dSNP/04.01.dSNP.RData")
dSNP.gene.drug.limma.mat = gene.drug.limma.mat; rm(gene.drug.limma.mat)
colnames(dSNP.gene.drug.limma.mat) = c("Cancer", "Gene", "Drug", "p.twosided", "FC", "nMut", "nMutProp")

dSNP.gene.drug.limma.mat = dSNP.gene.drug.limma.mat[!dSNP.gene.drug.limma.mat[,2] %in% long.genes, ]

dSNP.gene.drug.limma.mat[which(dSNP.gene.drug.limma.mat[,3]=="X17.AAG"),3] = "17-AAG"
gsub("\\.", "", dSNP.gene.drug.limma.mat[,3]) -> new.drug
dSNP.gene.drug.limma.mat[,3] = new.drug

dim(dSNP.gene.drug.limma.mat)

dSNP.gene.drug.limma.mat = dSNP.gene.drug.limma.mat[dSNP.gene.drug.limma.mat[,1] != "LGG" & dSNP.gene.drug.limma.mat[,1]!="BRCA" & dSNP.gene.drug.limma.mat[,1]!="THCA", ]
dim(dSNP.gene.drug.limma.mat)

write.table(dSNP.gene.drug.limma.mat, file="dSNP.gene.drug.txt", row.names=F, quote=F, sep="\t")
dSNP.gene.drug.limma.mat = read.table("dSNP.gene.drug.txt", as.is=T, header=T)

dSNP.gene.drug.limma.mat$adjp = NULL

###### *****
drugs = unique(dSNP.gene.drug.limma.mat[,3])
for(k in 1:length(drugs)){
	which(dSNP.gene.drug.limma.mat[,3] == drugs[k]) -> ii
	p.adjust(dSNP.gene.drug.limma.mat[ii,4], method="BH") -> adjp
	dSNP.gene.drug.limma.mat[ii, "adjp"] = adjp
}
###### *****


which(dSNP.gene.drug.limma.mat$adjp < 0.2 & dSNP.gene.drug.limma.mat$FC > 0) -> ii.1
pos.mat = dSNP.gene.drug.limma.mat[ii.1,]
which(dSNP.gene.drug.limma.mat$adjp < 0.2 & dSNP.gene.drug.limma.mat$FC < 0) -> ii.2
neg.mat = dSNP.gene.drug.limma.mat[ii.2,]
cc = rep(0, nrow(dSNP.gene.drug.limma.mat)); cc[ii.1] = 2; cc[ii.2] = 1
dSNP.gene.drug.limma.mat$color = as.factor(cc)
write.table(dSNP.gene.drug.limma.mat[which(cc!=0), ], file="dSNP.gene.drug.sig.txt", row.names=F, quote=F, sep="\t")
table(cc)


################################ !!!!!!!!!!!!!!!!!!!!!
library(ggrepel)
dat.plot = data.frame(log2FC = (as.numeric(dSNP.gene.drug.limma.mat[,5])), log10p=-log10(dSNP.gene.drug.limma.mat[,4]), nMutProp = as.numeric(dSNP.gene.drug.limma.mat[,7]), color = dSNP.gene.drug.limma.mat$color )

lab = apply(dSNP.gene.drug.limma.mat, 1, function(u)paste(u[1], u[2], u[3], sep="-") )
dat.plot$lab = lab

p2 = ggplot(dat.plot, aes(x=log2FC, y=log10p, size=nMutProp, color=color)) + geom_point(alpha=0.2) + scale_color_manual(values=c("black", "blue", "red")) + 
     theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size=8), axis.title.x = element_text(size=8) ) +
     geom_vline(xintercept = 0, linetype="dashed", color = "red", size=1) + ylab("-log10(p)") + xlab("Difference in means of MT and WT samples")
	 
which(dat.plot$log10p > 7.5) -> ii
new.dat.plot = dat.plot[ii, ]
p3 = p2+geom_text_repel(data=new.dat.plot, aes(label=lab), color="black", size=2)

tiff("6A.subtype.tif", width=1800, height=2000, res=300)
print(p3)
dev.off()
################################

##############################################################################################################

### plot3
which(-log10(dSNP.gene.drug.limma.mat$p.twosided) > 5 & -log10(dSNP.gene.drug.limma.mat$p.twosided) < 7.5 ) -> ii

dat.plot = data.frame(log2FC = (as.numeric(dSNP.gene.drug.limma.mat[ii,5])), log10p=-log10(dSNP.gene.drug.limma.mat[ii,4]), nMutProp = as.numeric(dSNP.gene.drug.limma.mat[ii,7]) )
which(dat.plot$log2FC < 0) -> ii.1
which(dat.plot$log2FC > 0) -> ii.2
cc = rep(0, nrow(dat.plot)); cc[ii.1] = 1; cc[ii.2] = 2
dat.plot$color = as.factor(cc)
levels(dat.plot$color) = c(0,1,2)

labels = apply(dSNP.gene.drug.limma.mat[ii,], 1, function(u)paste(u[1], u[2], u[3], sep="-") )
dat.plot$lab = labels

p3 = ggplot(dat.plot, aes(x=log2FC, y=log10p, size=nMutProp, color=color, label=lab)) + geom_point(alpha=0.2) + scale_color_manual(values=c("blue", "red")) + 
     theme(plot.title = element_text(hjust = 0.5, size=8), legend.position = c(0.1, 0.9), axis.title.x = element_text(size=8), legend.text=element_text(size=6), legend.title=element_text(size=6) ) + guides(color=F) +
	 geom_vline(xintercept = 0, linetype="dashed", color = "red", size=1) + ylab("-log10(p)") + xlab("Difference in means of MT and WT samples")

p3.2 = p3+geom_text_repel(data=dat.plot, aes(label=lab), color="black", size=1.5)


pdf("6B.subtype.pdf", width=5.5, height=5)
print(p3.2)
dev.off()

