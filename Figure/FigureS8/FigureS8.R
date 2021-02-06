#setwd("/path/to/VAEN/Figure/FigureS8")

load("../Figure3/3EF.data.RData")
sapply(rownames(TCGA.sensitive.mat), function(u)gsub("\\.", "-", u)) -> new.name
rownames(TCGA.sensitive.mat) = unlist(new.name)
sapply(rownames(TCGA.resistant.mat), function(u)gsub("\\.", "-", u)) -> new.name
rownames(TCGA.resistant.mat) = unlist(new.name)


TCGA.cutoff = 0.05

ii = 1:84
new.TCGA.sensitive.mat = TCGA.sensitive.mat[ii,]
new.TCGA.resistant.mat = TCGA.resistant.mat[ii,]

new.TCGA.sensitive.mat = new.TCGA.sensitive.mat[nrow(new.TCGA.sensitive.mat):1, ]
new.TCGA.resistant.mat = new.TCGA.resistant.mat[nrow(new.TCGA.resistant.mat):1, ]

pdf("S8.1.pdf", width=6, height=11)
par(mar=c(1,7,10,3))

m = -log10(new.TCGA.sensitive.mat)
m[which( abs(new.TCGA.sensitive.mat) > TCGA.cutoff )] = 0
m[which(is.na(m))] = 0
m = t(m)
m[which(m>100)] = 100

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2red <- colorRampPalette(c("white", "red"))
cc = c(white2red(500))


image(new.m, xaxt="n", yaxt="n", col = cc)

### row labels: drug names
geneCount = ncol(m)
path.col = rep("black", ncol(m))
lab = colnames(m)
lab[grep("Nutlin", lab)] = "Nutlin-3a (-)"
mtext(text=lab, side=2, line=0.3, at=(0:(geneCount-1))/(geneCount-1), las=1, cex=.8, col=path.col)

### col labels: cancer type
geneCount = nrow(m)
text((0:(geneCount-1))/(geneCount-1), 1.05, srt = 90, adj = 0, cex=.8, labels = rownames(m), xpd = TRUE)

geneCount = nrow(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(v=x, col=grey(0.9), lwd=.2)}
geneCount = ncol(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(h=x, col=grey(0.8), lwd=.2);}
box()

##############


m = -log10(new.TCGA.resistant.mat)
m[which( abs(new.TCGA.resistant.mat) > TCGA.cutoff )] = 0
m[which(is.na(m))] = 0
m = t(m)
m[which(m>100)] = 100

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2blue <- colorRampPalette(c("white", "blue"))
cc = c(white2blue(500))

geneCount = nrow(new.m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.x
geneCount = ncol(m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.y
for(k1 in 1:(length(tt.x)-1) ){
	for(k2 in 1:(length(tt.y)-1) ){
		rownames(m)[k1] -> gene
		colnames(m)[k2] -> pathway
		if(!(is.element(gene, rownames(new.m)) & is.element(pathway, colnames(new.m))))next
		if(new.m[gene, pathway]!=1){
			cc.idx = round(new.m[gene, pathway]/2 * 1e3)
			if(cc.idx < 1) cc.idx = 1
			polygon( c(tt.x[k1], tt.x[k1+1], tt.x[k1+1]), c(tt.y[k2], tt.y[k2], tt.y[k2+1]), col=cc[cc.idx], border=NA)
		}
	}
}
box()


mtext("(A)", side=2, at = 1.1, las=1, line=2)
dev.off()

#############################################################################

ii = 85:168
new.TCGA.sensitive.mat = TCGA.sensitive.mat[ii,]
new.TCGA.resistant.mat = TCGA.resistant.mat[ii,]

new.TCGA.sensitive.mat = new.TCGA.sensitive.mat[nrow(new.TCGA.sensitive.mat):1, ]
new.TCGA.resistant.mat = new.TCGA.resistant.mat[nrow(new.TCGA.resistant.mat):1, ]

pdf("S8.2.pdf", width=6, height=11)
par(mar=c(1,7,10,3))

m = -log10(new.TCGA.sensitive.mat)
m[which( abs(new.TCGA.sensitive.mat) > TCGA.cutoff )] = 0
m[which(is.na(m))] = 0
m = t(m)
m[which(m>100)] = 100

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2red <- colorRampPalette(c("white", "red"))
cc = c(white2red(500))


image(new.m, xaxt="n", yaxt="n", col = cc)

### row labels: drug names
geneCount = ncol(m)
path.col = rep("black", ncol(m))
lab = colnames(m)
lab[grep("Nutlin", lab)] = "Nutlin-3a (-)"
mtext(text=lab, side=2, line=0.3, at=(0:(geneCount-1))/(geneCount-1), las=1, cex=.8, col=path.col)

### col labels: cancer type
geneCount = nrow(m)
text((0:(geneCount-1))/(geneCount-1), 1.05, srt = 90, adj = 0, cex=.8, labels = rownames(m), xpd = TRUE)

geneCount = nrow(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(v=x, col=grey(0.9), lwd=.2)}
geneCount = ncol(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(h=x, col=grey(0.8), lwd=.2);}
box()

##############


m = -log10(new.TCGA.resistant.mat)
m[which( abs(new.TCGA.resistant.mat) > TCGA.cutoff )] = 0
m[which(is.na(m))] = 0
m = t(m)
m[which(m>100)] = 100

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2blue <- colorRampPalette(c("white", "blue"))
cc = c(white2blue(500))

geneCount = nrow(new.m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.x
geneCount = ncol(m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.y
for(k1 in 1:(length(tt.x)-1) ){
	for(k2 in 1:(length(tt.y)-1) ){
		rownames(m)[k1] -> gene
		colnames(m)[k2] -> pathway
		if(!(is.element(gene, rownames(new.m)) & is.element(pathway, colnames(new.m))))next
		if(new.m[gene, pathway]!=1){
			cc.idx = round(new.m[gene, pathway]/2 * 1e3)
			if(cc.idx < 1) cc.idx = 1
			polygon( c(tt.x[k1], tt.x[k1+1], tt.x[k1+1]), c(tt.y[k2], tt.y[k2], tt.y[k2+1]), col=cc[cc.idx], border=NA)
		}
	}
}
box()

mtext("(B)", side=2, at = 1.1, las=1, line=2)
dev.off()



#############################################################################

ii = 169:251
new.TCGA.sensitive.mat = TCGA.sensitive.mat[ii,]
new.TCGA.resistant.mat = TCGA.resistant.mat[ii,]

new.TCGA.sensitive.mat = new.TCGA.sensitive.mat[nrow(new.TCGA.sensitive.mat):1, ]
new.TCGA.resistant.mat = new.TCGA.resistant.mat[nrow(new.TCGA.resistant.mat):1, ]

pdf("S8.3.pdf", width=6, height=11)
par(mar=c(1,7,10,3))

m = -log10(new.TCGA.sensitive.mat)
m[which( abs(new.TCGA.sensitive.mat) > TCGA.cutoff )] = 0
m[which(is.na(m))] = 0
m = t(m)
m[which(m>100)] = 100

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2red <- colorRampPalette(c("white", "red"))
cc = c(white2red(500))


image(new.m, xaxt="n", yaxt="n", col = cc)

### row labels: drug names
geneCount = ncol(m)
path.col = rep("black", ncol(m))
lab = colnames(m)
lab[grep("Nutlin", lab)] = "Nutlin-3a (-)"
mtext(text=lab, side=2, line=0.3, at=(0:(geneCount-1))/(geneCount-1), las=1, cex=.8, col=path.col)

### col labels: cancer type
geneCount = nrow(m)
text((0:(geneCount-1))/(geneCount-1), 1.05, srt = 90, adj = 0, cex=.8, labels = rownames(m), xpd = TRUE)

geneCount = nrow(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(v=x, col=grey(0.9), lwd=.2)}
geneCount = ncol(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(h=x, col=grey(0.8), lwd=.2);}
box()


m = -log10(new.TCGA.resistant.mat)
m[which( abs(new.TCGA.resistant.mat) > TCGA.cutoff )] = 0
m[which(is.na(m))] = 0
m = t(m)
m[which(m>100)] = 100

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2blue <- colorRampPalette(c("white", "blue"))
cc = c(white2blue(500))

geneCount = nrow(new.m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.x
geneCount = ncol(m)
-1/(2*geneCount-2) + (0:(geneCount))/(geneCount-1) -> tt.y
for(k1 in 1:(length(tt.x)-1) ){
	for(k2 in 1:(length(tt.y)-1) ){
		rownames(m)[k1] -> gene
		colnames(m)[k2] -> pathway
		if(!(is.element(gene, rownames(new.m)) & is.element(pathway, colnames(new.m))))next
		if(new.m[gene, pathway]!=1){
			cc.idx = round(new.m[gene, pathway]/2 * 1e3)
			if(cc.idx < 1) cc.idx = 1
			polygon( c(tt.x[k1], tt.x[k1+1], tt.x[k1+1]), c(tt.y[k2], tt.y[k2], tt.y[k2+1]), col=cc[cc.idx], border=NA)
		}
	}
}
box()

mtext("(C)", side=2, at = 1.1, las=1, line=2)


############## legend
##############

col = cc
breaks = seq(0,500,1)/500
poly <- vector(mode="list", length(col))
for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
}
  
for(i in seq(poly)){
  polygon(c(0.9,0.9,0.95,0.95), poly[[i]], col=col[i], border=NA)
}

##############
##############
cc = c(white2red(500))

col = cc
breaks = seq(0,500,1)/500
poly <- vector(mode="list", length(col))
for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
}
  
for(i in seq(poly)){
  polygon(c(0.95,0.95,1,1), poly[[i]], col=col[i], border=NA)
}

##############
##############


dev.off()
