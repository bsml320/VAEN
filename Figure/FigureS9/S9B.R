setwd("/path/to/VAEN/Figure/FigureS9")
library("reshape2")
library("ggplot2")
library("plotrix")

############################################################

load("../Figure3/3EF.data.RData")

blue2white <- colorRampPalette(c("blue", "white"))
white2red <- colorRampPalette(c("white", "red"))
cc = c(blue2white(100), white2red(100))

##################################################################
##################################################################
two.drugs.match = read.delim("../../DATA/drugs.match-2.txt", as.is=T, header=F)
two.drugs.match = two.drugs.match[order(two.drugs.match[,3]), ]

match(two.drugs.match[,3], rownames(obsd.drug.by.tissue.mat)) -> shared.ii
shared.ii = shared.ii[length(shared.ii):1]

which( apply(obsd.drug.by.tissue.mat, 1, function(u)sum(abs(u) < 0.05)) > 8) -> ii
unique.ii = setdiff(ii, shared.ii)
unique.ii = unique.ii[length(unique.ii):1]
draw.ii = c( unique.ii, shared.ii )


obsd.mat = obsd.drug.by.tissue.mat[draw.ii, ]
pred.mat = pred.drug.by.tissue.mat[draw.ii, ]
full.mat = full.drug.by.tissue.mat[draw.ii, ]

x = -log10(abs(obsd.mat)) * sign(obsd.mat)
y = -log10(abs(pred.mat)) * sign(pred.mat)
z = -log10(abs(full.mat)) * sign(full.mat)

pos = max(  c( max(x[x>0]), max(y[y>0]), max(z[z>0]) ) )
neg = min(  c( min(x[x<0]), min(y[y<0]), min(z[z<0]) ) )

############################################################

m = -log10(abs(obsd.mat)) * sign(obsd.mat)
m[which( abs(obsd.mat) > 0.05 )] = 0
m[which(is.na(m))] = 0
m = t(m)

m1 = m; m1[which(m<0)] = 0;
m2 = m; m2[which(m>0)] = 0; m2 = abs(m2)

m11 = m1/pos
m22 = m2/abs(neg)

obsd.m = m11 + (-m22)

############################################################

rppa.m = -log10(abs(pred.mat)) * sign(pred.mat)
rppa.m[which( abs(pred.mat) > 0.05 )] = 0
rppa.m[which(is.na(rppa.m))] = 0
rppa.m = t(rppa.m)

m1 = rppa.m; m1[which(rppa.m<0)] = 0;
m2 = rppa.m; m2[which(rppa.m>0)] = 0; m2 = abs(m2)

m11 = m1/pos
m22 = m2/abs(neg)

pred.m = m11 + (-m22)

############################################################

rppa.m = -log10(abs(full.mat)) * sign(full.mat)
rppa.m[which( abs(full.mat) > 0.05 )] = 0
rppa.m[which(is.na(rppa.m))] = 0
rppa.m = t(rppa.m)

m1 = rppa.m; m1[which(rppa.m<0)] = 0;
m2 = rppa.m; m2[which(rppa.m>0)] = 0; m2 = abs(m2)

m11 = m1/pos
m22 = m2/abs(neg)

full.m = m11 + (-m22)

############## legend
pdf("FigureS9B.legend.pdf")

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

########################################################################################################

pdf("FigureS9B.pdf", width=4.2, height=4.8)
par(mar=c(0.5,6.8,6,4))

image(obsd.m, xaxt="n", yaxt="n", col = "white")

geneCount = ncol(m)
path.col = rep("black", ncol(m))
colnames(obsd.m) -> label
label[which(label == "Crizotinib")] = "PF2341066 | Crizotinib"
label[which(label == "NVP-TAE684")] = "TAE684 | NVP-TAE684"
label[which(label == "Palbociclib")] = "PD-0332991 | Palbociclib"
label[which(label == "Saracatinib")] = "AZD0530 | Saracatinib"
label[which(label == "Tanespimycin")] = "17-AAG | Tanespimycin"
mtext(text=label, side=2, line=0.3, at=(0:(geneCount-1))/(geneCount-1), las=1, cex=.7, col=path.col)


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
mtext("(C)", side=2, at = 1.1, las=1, line=2, cex=.9)
mtext("(D)", side=2, at = 1.2, las=1, line=2, cex=.9)

dev.off()
