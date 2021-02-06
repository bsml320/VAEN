setwd("/path/to/VAEN/Figure/FigureS9")

load("../Figure3/3BCD.data.RData")

blue2white <- colorRampPalette(c("blue", "white"))
white2red <- colorRampPalette(c("white", "red"))
cc = c(blue2white(100), white2red(100))

obsd.drug.by.tissue.mat = obsd.drug.by.tissue.mat[nrow(obsd.drug.by.tissue.mat):1, ]
pred.drug.by.tissue.mat = pred.drug.by.tissue.mat[nrow(pred.drug.by.tissue.mat):1, ]
full.drug.by.tissue.mat = full.drug.by.tissue.mat[nrow(full.drug.by.tissue.mat):1, ]

##############
x = -log10(abs(full.drug.by.tissue.mat)) * sign(full.drug.by.tissue.mat)
y = -log10(abs(obsd.drug.by.tissue.mat)) * sign(obsd.drug.by.tissue.mat)
z = -log10(abs(pred.drug.by.tissue.mat)) * sign(pred.drug.by.tissue.mat)

pos = max(  c( max(x[x>0]), max(y[y>0]), max(z[z>0]) ) )
neg = min(  c( min(x[x<0]), min(y[y<0]), min(z[z<0]) ) )

############## legend
pdf("FigureS9A.legend.pdf")

plot(1,1, xlim=c(0,1), ylim=c(0,1))
col = white2red(100)
breaks = seq(0,100,1)/100
poly <- vector(mode="list", length(col))
for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
}
  
for(i in seq(poly)){
  polygon(poly[[i]], c(0.1,0.1,0.15,0.15), col=col[i], border=NA)
}

text(0.1,0.1,"-log10(p), 0")
text(poly[[i]][1], 0.1 ,round(pos))


plot(-1,0, xlim=c(-1,0), ylim=c(0,1))
col = blue2white(100)
breaks = seq(-100,0,1)/100
poly <- vector(mode="list", length(col))
for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
}
  
for(i in seq(poly)){
  polygon(poly[[i]], c(-1,-1,0,0), col=col[i], border=NA)
}

text(0.1,0.1,"-log10(p), 0")
text(-1, 0.1 ,round(neg))


dev.off()

######################################################################

m = -log10(abs(full.drug.by.tissue.mat)) * sign(full.drug.by.tissue.mat)
m[which( abs(full.drug.by.tissue.mat) > 0.05 )] = 0
m[which(is.na(m))] = 0
m = t(m)

m1 = m; m1[which(m<0)] = 0;
m2 = m; m2[which(m>0)] = 0; m2 = abs(m2)
m11 = m1/pos
m22 = m2/abs(neg)
full.m = m11 + (-m22)

##############

rppa.m = -log10(abs(pred.drug.by.tissue.mat)) * sign(pred.drug.by.tissue.mat)
rppa.m[which( abs(pred.drug.by.tissue.mat) > 0.05 )] = 0
rppa.m[which(is.na(rppa.m))] = 0
rppa.m = t(rppa.m)

m1 = rppa.m; m1[which(rppa.m<0)] = 0;
m2 = rppa.m; m2[which(rppa.m>0)] = 0; m2 = abs(m2)
m11 = m1/pos
m22 = m2/abs(neg)
pred.m = m11 + (-m22)

##############

rppa.m = -log10(abs(obsd.drug.by.tissue.mat)) * sign(obsd.drug.by.tissue.mat)
rppa.m[which( abs(obsd.drug.by.tissue.mat) > 0.05 )] = 0
rppa.m[which(is.na(rppa.m))] = 0
rppa.m = t(rppa.m)

m1 = rppa.m; m1[which(rppa.m<0)] = 0;
m2 = rppa.m; m2[which(rppa.m>0)] = 0; m2 = abs(m2)
m11 = m1/pos
m22 = m2/abs(neg)
obsd.m = m11 + (-m22)

##########################################################################################################

library("plotrix")

blue2white <- colorRampPalette(c("blue", "white"))
white2red <- colorRampPalette(c("white", "red"))
cc = c(blue2white(100), white2red(100))

pdf("FigureS9A.pdf", width=4.18, height=5)
par(mar=c(2,6.8,10,4))

image(obsd.m, xaxt="n", yaxt="n", col = "white")
geneCount = ncol(m)
path.col = rep("black", ncol(m))
mtext(text=colnames(m), side=2, line=0.3, at=(0:(geneCount-1))/(geneCount-1), las=1, cex=.7, col=path.col)
geneCount = nrow(m)

text((0:(geneCount-1))/(geneCount-1), 1.05, srt = 45, adj = 0, cex=.5, labels = rownames(m), xpd = TRUE)

geneCount = nrow(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(v=x, col=grey(0.9), lwd=.2)}
geneCount = ncol(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(h=x, col=grey(0.8), lwd=.2);}


geneCount = nrow(m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.x
row.length = (tt.x[2] - tt.x[1])
geneCount = ncol(m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.y
for(k1 in 1:(length(tt.x)-1) ){
	for(k2 in 1:(length(tt.y)-1) ){
		rownames(m)[k1] -> gene
		colnames(m)[k2] -> pathway
		if(!(is.element(gene, rownames(rppa.m)) & is.element(pathway, colnames(rppa.m))))next
		if(rppa.m[gene, pathway]!=2){
			obsd.idx = 100 + round(obsd.m[gene, pathway] * 1e2)
			pred.idx = 100 + round(pred.m[gene, pathway] * 1e2)
			full.idx = 100 + round(full.m[gene, pathway] * 1e2)
			
			xc = tt.x[k1] + row.length/2
			yc = (tt.y[k2] + tt.y[k2+1])/2
			
			floating.pie( xc, yc, c(1,1,1), col=c(cc[obsd.idx],cc[pred.idx],cc[full.idx]),radius=row.length/2 * sin(1), border=NA )  
		}
	}
}
box()

mtext("(A)", side=2, at = 1.1, las=1, line=2, cex=.9)
mtext("(B)", side=2, at = 1.2, las=1, line=2, cex=.9)

dev.off()

pdf("legend.pdf")
par(mar=rep(10, 4))
plot(1)
floating.pie( 1, 1, c(1,1,1), col=c("red","orange","pink"),radius=0.7/2 * sin(1), border=NA )
mtext("Set 1 (obsd)", side=2, at = 1.1, las=1, line=2, cex=.9)
mtext("Set 2 (pred)", side=2, at = 1.2, las=1, line=2, cex=.9)
mtext("Set 3 (full)", side=2, at = 1.3, las=1, line=2, cex=.9)
dev.off()
