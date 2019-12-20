setwd("/work/Figures/Figure6/")
source("/work/code/multiplot.R")

load("/work/data/UCSC/ncbiRefSeq.01252019.hg19.gene2ll.RData")
long.genes = names(gene2ll[gene2ll > 200000])

#################################################################################################

dSNP.gene.drug.limma.mat = read.table("dSNP.gene.drug.txt", as.is=T, header=T)  ### from 6AB.subtype.R
dSNP.sig = read.delim("dSNP.gene.drug.sig.txt", as.is=T)
dSNP.sig = dSNP.sig[!dSNP.sig[,2] %in% long.genes, ]
drugs = unique(dSNP.gene.drug.limma.mat[,3])

#####################
table(dSNP.sig[,2]) -> count
rec.genes = names(which(count > 2))

dSNP.sig = dSNP.sig[which(dSNP.sig[,2] %in% rec.genes), ]

apply(dSNP.sig, 1, function(u)paste(u[2], u[1], sep="-") ) -> gene.cancer.pair
gene.cancer.pair.unique = sort( unique(gene.cancer.pair) )

mut.mat = matrix(1, nrow=length(gene.cancer.pair.unique), ncol=length(drugs))
rownames(mut.mat) = gene.cancer.pair.unique
colnames(mut.mat) = drugs

for(k in 1:nrow(mut.mat)){
	x = strsplit(gene.cancer.pair.unique[k], split="-")[[1]]
	gene   = strsplit(gene.cancer.pair.unique[k], split="-")[[1]][1]
	cancer = strsplit(gene.cancer.pair.unique[k], split="-")[[1]][2]
	if(length(x) > 2){
		cancer = paste(strsplit(gene.cancer.pair.unique[k], split="-")[[1]][-1], collapse="-")
	}
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		which(dSNP.gene.drug.limma.mat[,1]==cancer & dSNP.gene.drug.limma.mat[,2] == gene & dSNP.gene.drug.limma.mat[,3]==drug) -> ii
		if(length(ii) > 0)mut.mat[k, kdrug] = dSNP.gene.drug.limma.mat[ii, 4] * sign(dSNP.gene.drug.limma.mat[ii, 5])
	}
}

##############################################
library(corrplot)

m = -log10(abs(mut.mat)) * sign(mut.mat)
m = m/max(abs(m))
col1 = colorRampPalette(c("blue", "white"))
col2 = colorRampPalette(c("white", "red"))
col = c(col1(90), rep("white", 20), col2(90))
col = c(col1(100), col2(100))

min(m[m!=0])
max(m)

pdf("S.6C.subtype.pdf")
corrplot(m, col=col)
dev.off()

##############################################
##############################################
new.mut.mat = apply(mut.mat, 2, function(u){
	sign(u) -> lab
	which(abs(u) < 1e-8) -> ii
	u[ii] = 1e-8 * lab[ii]
	u
})

m = -log10(abs(new.mut.mat)) * sign(new.mut.mat)
m = m/max(abs(m))
col1 = colorRampPalette(c("blue", "white"))
col2 = colorRampPalette(c("white", "red"))
col = c(col1(90), rep("white", 20), col2(90))
col = c(col1(100), col2(100))

min(m[m!=0])
max(m)

pdf("6C.scaled.subtype.pdf", width=6, height=8)
corrplot(m, col=col)
dev.off()



