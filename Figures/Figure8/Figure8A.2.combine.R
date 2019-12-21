setwd("C:/Users/pjia/UTH/road/18-VAE/analysis/V15/result.EN/dr.CCLE/04-mix/09.Expr/DAVID")
ff = list.files("David_enrichment_plot_GOBP")
ff = ff[grep("txt", ff)]

go.bp = c()
for(k in 1:length(ff)){
	dat = read.delim(paste("David_enrichment_plot_GOBP/", ff[k], sep=""), as.is=T)
	which(dat[,2] >= 10 & dat[,3] < 1e-8) -> ii
	if(length(ii) > 0)go.bp = c(go.bp, dat[ii,1])
}
go.bp = unique(go.bp)
length(go.bp)

bp.mat = matrix(1, nrow=length(go.bp), ncol=length(ff))
rownames(bp.mat) = go.bp
f1 = gsub(".Q95.txt", "", ff)
f1 = gsub(".high", ".H", f1)
f1 = gsub(".low", ".L", f1)
colnames(bp.mat) = f1

for(k in 1:length(ff)){
	dat = read.delim(paste("David_enrichment_plot_GOBP/", ff[k], sep=""), as.is=T)
	which(dat[,2] >= 10 & dat[,3] < 1e-8) -> ii
	if(length(ii) > 0){
		match(dat[ii, 1], go.bp) -> idx
		direction = 1
		if(grepl("low", ff[k]))direction = -1
		bp.mat[idx, k] = dat[ii, 3] * direction
	}
}

apply(bp.mat, 1, function(u)sum( abs(u) < 1e-8)) -> check
new.bp.mat = bp.mat[which(check >= 2), ]

heatmap.2(new.bp.mat, trace="none") -> hc
new.dat = new.bp.mat[hc$rowInd, hc$colInd]

slim.new.dat = new.dat
for(k1 in 1:nrow(new.dat)){
	for(k2 in 1:ncol(new.dat)){
		slim.new.dat[k1, k2] = sign(new.dat[k1, k2]) * (-log10( abs(new.dat[k1, k2]))  )
	}
}
new.dat = slim.new.dat
heatmap.2(new.dat, trace="none") -> hc

scale.factor = max(abs(new.dat))
new.dat[new.dat > 50] = 50
new.dat = new.dat/scale.factor

for(k in 1:length(ff)){
	dat = read.delim(paste("David_enrichment_plot_GOBP/", ff[k], sep=""), as.is=T)
	intersect(dat[,1], rownames(new.dat)) -> ii
	if(length(ii) > 0){
		match(intersect(dat[,1], rownames(new.dat)), rownames(new.dat)) -> idx
		direction = 1
		if(grepl("low", ff[k]))direction = -1
		new.dat[idx, k] = (-log10(dat[ match(ii, dat[,1]) , 3]) ) * direction / scale.factor
	}
}
corrplot(new.dat, method="ellipse")


pdf("DAVID.BP.pdf", width=14, height=10)
corrplot(new.dat, method="ellipse")
dev.off()

p.mat = 10^(-abs(new.dat) * scale.factor)
pdf("DAVID.BP.pdf", width=14, height=10)
corrplot(new.dat, p.mat = p.mat, method="ellipse", sig.level=1e-8, pch.cex = .75)
dev.off()

p.mat2 = p.mat
p.mat2[p.mat > 1e-8] <- 1e-8
p.mat2[p.mat < 1e-8] <- 1

pdf("DAVID.BP.pdf", width=14, height=10)
corrplot(new.dat, p.mat = p.mat2, method="ellipse", sig.level=1e-8, pch.cex = .75)
dev.off()
