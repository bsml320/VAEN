#setwd("/path/to/VAEN/Figure/Figure7")

library("ggplot2")
library("ggrepel")
library("corrplot")

#################################################################################################
#################################################################################################

dSNP.cluster.mat = read.table("Figure7C.txt", header=T, as.is=T)

print(dSNP.cluster.mat[dSNP.cluster.mat[,2]=="FLT3",])

sapply(as.character(dSNP.cluster.mat[,8]), function(u){strsplit(u, split=":")[[1]] -> v; 
	if(length(unique(v)) > 1) {
		paste(min(as.numeric(v)), max(as.numeric(v)), sep="-") -> v1;
	} else {
		unique(v) -> v1;
	}
	v1
}) -> poss
dSNP.cluster.mat$poss = poss


dSNP.sig = read.delim("Figure7C.sig.txt", as.is=T)
dSNP.sig = dSNP.sig[dSNP.sig[,2]!="TP53", ]   ### excluding TP53 from plotting

apply(dSNP.sig, 1, function(u)paste(u[1], u[2], u[4], u[5], collapse=":")) -> lab
match(unique(lab), lab) -> ii
dSNP.sig = dSNP.sig[ii, ]

#####################
table(dSNP.sig[,2]) -> count
rec.genes = names(which(count > 1))

dSNP.sig = dSNP.sig[which(dSNP.sig[,2] %in% rec.genes), ]

sapply(as.character(dSNP.sig[,8]), function(u){strsplit(u, split=":")[[1]] -> v; 
	if(length(unique(v)) > 1) {
		paste(min(as.numeric(v)), max(as.numeric(v)), sep="-") -> v1;
	} else {
		unique(v) -> v1;
	}
	v1
}) -> poss
dSNP.sig$poss = poss

apply(dSNP.sig, 1, function(u)paste(u[2], u[3], u[1], u[12], sep=":") ) -> gene.cancer.pair
gene.cancer.pair.unique = sort( unique(gene.cancer.pair) )

mut.mat = matrix(1, nrow=length(gene.cancer.pair.unique), ncol=length(drugs))
rownames(mut.mat) = gene.cancer.pair.unique
colnames(mut.mat) = drugs

for(k in 1:nrow(mut.mat)){
	gene   = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][1]
	trans  = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][2]
	cancer = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][3]
	poss   = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][4]
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		which(dSNP.cluster.mat[,1] == cancer & dSNP.cluster.mat[,2] == gene & dSNP.cluster.mat[,3] == trans & dSNP.cluster.mat[,12] == poss & dSNP.cluster.mat[,4]==drug) -> ii
		if(length(ii) == 0)next
		which.min(dSNP.cluster.mat[ii, 5]) -> idx
		mut.mat[k, kdrug] = dSNP.cluster.mat[ii[idx], 5] * ifelse(dSNP.cluster.mat[ii[idx], 6] > 0, 1, -1)
	}
}

sapply(rownames(mut.mat), function(u){
	strsplit(u, split=":")[[1]] -> v
	paste(v[3], v[1], v[4], sep=":")
}) -> new.name
rownames(mut.mat) = new.name

##############################################

m = -log10(abs(mut.mat)) * sign(mut.mat)

min(m[m!=0])
max(m)

x = ceiling(max(abs(m)))

pdf("7C.pdf", height=12)
corrplot(m,is.corr = F,cl.align.text = "l",cl.offset = 0.3,cl.lim = c(-x,x),col=colorRampPalette(c("blue","white","red"))(256))
dev.off()

