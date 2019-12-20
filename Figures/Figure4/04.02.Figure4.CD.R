setwd("/work/Figures/Figure4/")
library("reshape2")
load("/work/Figures/Figure3/GDSC.lineage.RData")   ### use the data generated in Figure 3

##################################################################
two.drugs.match = read.delim("/work/data/drugs.match-2.txt", as.is=T, header=F)
two.drugs.match = two.drugs.match[order(two.drugs.match[,3]), ]

match(two.drugs.match[,3], rownames(obsd.drug.by.tissue.mat)) -> shared.ii
shared.ii = shared.ii[length(shared.ii):1]


##################################################################
############ Figure 4C

library("plotrix")

blue2white <- colorRampPalette(c("blue", "white"))
white2red <- colorRampPalette(c("white", "red"))
cc = c(blue2white(100), white2red(100))

### prepare matrix for observed drug response in GDSC
which( apply(obsd.drug.by.tissue.mat, 1, function(u)sum(abs(u) < 0.05)) > 8) -> ii
unique.ii = setdiff(ii, shared.ii)
unique.ii = unique.ii[length(unique.ii):1]
draw.ii = c( unique.ii, shared.ii )
obsd.mat = obsd.drug.by.tissue.mat[draw.ii, ]

m = -log10(abs(obsd.mat)) * sign(obsd.mat)
m[which( abs(obsd.mat) > 0.05 )] = 0
m[which(is.na(m))] = 0
m = t(m)

m1 = m; m1[which(m<0)] = 0;
m2 = m; m2[which(m>0)] = 0; m2 = abs(m2)
m11 = (m1-min(m1, na.rm=T))/(max(m1, na.rm=T) - min(m1, na.rm=T))
m22 = (m2-min(m2, na.rm=T))/(max(m2, na.rm=T) - min(m2, na.rm=T))
obsd.m = m11 + (-m22)

### prepare matrix for predicted drug response by GDSC models
pred.mat = pred.drug.by.tissue.mat[draw.ii, ]
raw.m = -log10(abs(pred.mat)) * sign(pred.mat)
raw.m[which( abs(pred.mat) > 0.05 )] = 0
raw.m[which(is.na(raw.m))] = 0
raw.m = t(raw.m)

m1 = raw.m; m1[which(raw.m<0)] = 0;
m2 = raw.m; m2[which(raw.m>0)] = 0; m2 = abs(m2)
m11 = (m1-min(m1, na.rm=T))/(max(m1, na.rm=T) - min(m1, na.rm=T))
m22 = (m2-min(m2, na.rm=T))/(max(m2, na.rm=T) - min(m2, na.rm=T))
pred.m = m11 + (-m22)

### prepare matrix for predicted drug response by GDSC models in the full set of cell lines
full.mat = full.drug.by.tissue.mat[draw.ii, ]
rppa.m = -log10(abs(full.mat)) * sign(full.mat)
rppa.m[which( abs(full.mat) > 0.05 )] = 0
rppa.m[which(is.na(rppa.m))] = 0
rppa.m = t(rppa.m)

m1 = rppa.m; m1[which(rppa.m<0)] = 0;
m2 = rppa.m; m2[which(rppa.m>0)] = 0; m2 = abs(m2)
m11 = (m1-min(m1, na.rm=T))/(max(m1, na.rm=T) - min(m1, na.rm=T))
m22 = (m2-min(m2, na.rm=T))/(max(m2, na.rm=T) - min(m2, na.rm=T))
full.m = m11 + (-m22)


pdf("Figure4C.lineage.pdf", width=4.2, height=4.8)
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
		if(!(is.element(gene, rownames(raw.m)) & is.element(pathway, colnames(raw.m))))next
		if(raw.m[gene, pathway]!=2){
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


##################################################################
############ Figure 4D

rownames(TCGA.drug.by.cancer.mat) = rownames(obsd.drug.by.tissue.mat)

draw.mat = TCGA.drug.by.cancer.mat[draw.ii,]
rownames(draw.mat) = label
melt(draw.mat) -> mm
x = mm[,3]

###
which( abs(mm[,3]) < 1e-100) -> ii
mm[ii,3] = sign(x)[ii] * 1e-100
x = mm[, 3]
###

mm[,3] = -log10(abs(x)) * sign(x)

ggheatmap.3 <- ggplot(data = mm, aes(Var2, Var1, fill = value))+
 geom_tile(color = "grey")+
 scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 0, space = "Lab", name="t-test") +
 theme_minimal()+ 
 theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 0), axis.text.y = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank() )+ scale_x_discrete(position="top")+
 coord_fixed()

pdf("Figure4D.lineage.pdf", width=6.2, height=4.1)
print(ggheatmap.3)
dev.off()

