setwd("/work/Figures/Figure6/")
library("ggplot2")
library("ggrepel")

dSNP.cluster.mat = read.table("/work/result.EN/dr.CCLE/04-mix/05.cluster.dSNP/04.05.cluster.dSNP.txt", header=F, as.is=T, sep="\t")
dSNP.cluster.mat = dSNP.cluster.mat[, -5]
dSNP.cluster.mat = dSNP.cluster.mat[dSNP.cluster.mat[,1] != "LGG" & dSNP.cluster.mat[,1]!="BRCA" & dSNP.cluster.mat[,1]!="THCA", ]

grep(";", dSNP.cluster.mat[,2]) -> ii
dSNP.cluster.mat = dSNP.cluster.mat[-ii, ]

################### update drug name
which(dSNP.cluster.mat[,4] == "X17.AAG") -> ii
dSNP.cluster.mat[ii,4] = "17-AAG"
gsub("\\.", "", dSNP.cluster.mat[,4]) -> ss
dSNP.cluster.mat[,4] = ss 
###############################################################################################

dSNP.cluster.ps = as.numeric(dSNP.cluster.mat[,5])

###### *****
drugs = unique(dSNP.cluster.mat[,4])
dSNP.cluster.mat$adjp = 1
for(k in 1:length(drugs)){
	which(dSNP.cluster.mat[,4] == drugs[k]) -> ii
	p.adjust(dSNP.cluster.mat[ii,5], method="BH") -> adjp
	dSNP.cluster.mat[ii, "adjp"] = adjp
}
###### *****


which(dSNP.cluster.mat$adjp < 0.2 & dSNP.cluster.mat[,6] > 0) -> ii.1
pos.mat = dSNP.cluster.mat[ii.1,]
which(dSNP.cluster.mat$adjp < 0.2 & dSNP.cluster.mat[,6] < 0) -> ii.2
neg.mat = dSNP.cluster.mat[ii.2,]
cc = rep(0, nrow(dSNP.cluster.mat)); cc[ii.1] = 2; cc[ii.2] = 1
dSNP.cluster.mat$color = as.factor(cc)
table(cc)
sum(cc!=0)
write.table(dSNP.cluster.mat[which(cc!=0), ], file="dSNP.cluster.sig.txt", row.names=F, quote=F, sep="\t")
write.table(dSNP.cluster.mat, file="dSNP.cluster.txt", row.names=F, quote=F, sep="\t")

#################################################################################################
dSNP.cluster.mat = read.table("dSNP.cluster.txt", header=T, as.is=T)
tapply(dSNP.cluster.mat[,3], dSNP.cluster.mat[,2], function(u){names(which.max(table(u)))})  -> x

dSNP.sig = read.delim("dSNP.cluster.sig.txt", as.is=T)
dSNP.sig = dSNP.sig[!dSNP.sig[,2] %in% long.genes, ]

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

apply(dSNP.sig, 1, function(u)paste(u[2], u[1], u[12], sep=":") ) -> gene.cancer.pair
gene.cancer.pair.unique = sort( unique(gene.cancer.pair) )

mut.mat = matrix(1, nrow=length(gene.cancer.pair.unique), ncol=length(drugs))
rownames(mut.mat) = gene.cancer.pair.unique
colnames(mut.mat) = drugs

for(k in 1:nrow(mut.mat)){
	gene   = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][1]
	cancer = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][2]
	poss   = strsplit(gene.cancer.pair.unique[k], split=":")[[1]][3]
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		which(dSNP.sig[,1] == cancer & dSNP.sig[,2] == gene & dSNP.sig[,12] == poss & dSNP.sig[,4]==drug) -> ii
		if(length(ii) == 0)next
		which.min(dSNP.sig[ii, 5]) -> idx
		mut.mat[k, kdrug] = dSNP.sig[ii[idx], 5] * ifelse(dSNP.sig[ii[idx], 6] > 0, 1, -1)
	}
}

##############################################

m = -log10(abs(mut.mat)) * sign(mut.mat)
m = m/max(abs(m))
col2 = colorRampPalette(c("blue", "white", "red"))

min(m[m!=0])
max(m)

pdf("dSNP.cluster.sig.all.pdf", height=20)
col = col2(100)
corrplot(m, col=col)
dev.off()


#################################################################################################
dSNP.cluster.mat = read.table("dSNP.cluster.txt", header=T, as.is=T)

sapply(as.character(dSNP.cluster.mat[,8]), function(u){strsplit(u, split=":")[[1]] -> v; 
	if(length(unique(v)) > 1) {
		paste(min(as.numeric(v)), max(as.numeric(v)), sep="-") -> v1;
	} else {
		unique(v) -> v1;
	}
	v1
}) -> poss
dSNP.cluster.mat$poss = poss


dSNP.sig = read.delim("dSNP.cluster.sig.txt", as.is=T)
dSNP.sig = dSNP.sig[!dSNP.sig[,2] %in% long.genes, ]

apply(dSNP.sig, 1, function(u)paste(u[1], u[2], u[4], u[5], collapse=":")) -> lab
match(unique(lab), lab) -> ii
dSNP.sig = dSNP.sig[ii, ]

dSNP.sig = dSNP.sig[dSNP.sig[,2]!="TP53", ]

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
p.mat = abs(mut.mat)
m = -log10(abs(mut.mat)) * sign(mut.mat)
m = m/max(abs(m))
col2 = colorRampPalette(c("blue", "white", "red"))
col1 = colorRampPalette(c("blue", "white"))
col2 = colorRampPalette(c("white", "red"))
col = c( col1(110)[1:95], rep("white", 10), col2( 110 )[-c(1:15)])

min(m[m!=0])
max(m)

pdf("6D.dSNP.cluster.sig.no.TP53.pdf", height=8)
corrplot(m, col=col)
dev.off()


