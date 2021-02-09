#setwd("/path/to/VAEN/Figure/Figure7")

#################################################################################################

dSNP.gene.drug.limma.mat = read.table("Figure7A.txt", as.is=T, header=T)
dSNP.sig = read.delim("Figure7B.txt", as.is=T)

drugs = unique(dSNP.gene.drug.limma.mat[,3])

#####################
table(dSNP.sig[,2]) -> count
rec.genes = names(which(count > 2))
rec.genes = setdiff(rec.genes, "TP53")

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
x = ceiling(max(abs(m)))

pdf("7B.pdf", width=8, height=10)
corrplot(m,is.corr = F,cl.align.text = "l",cl.offset = 0.3,cl.lim = c(-x,x),col=colorRampPalette(c("blue","white","red"))(256))
dev.off()


