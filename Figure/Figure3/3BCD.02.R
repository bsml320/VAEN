#setwd("/path/to/VAEN/Figure/Figure3")
load("3BCD.data.RData")

pdf("3C.CCLE.lineage.pdf", width=5.7, height=5)
TCGA.cutoff = 0.05
par(mar=c(1,6,10,4))

TCGA.sensitive.mat = TCGA.sensitive.mat[nrow(TCGA.sensitive.mat):1, ]
TCGA.resistant.mat = TCGA.resistant.mat[nrow(TCGA.resistant.mat):1, ]

m = -log10(TCGA.sensitive.mat)

which(is.infinite(m)) -> ii
if(length(ii)>0){
	m[which(is.infinite(m))] = max(m[!is.infinite(m)])
}

m[which( abs(m) < -log10(TCGA.cutoff) )] = 0
m[which( abs(m) > 100 )] = 100
m[which(is.na(m))] = 0
m = t(m)

new.m = (m-min(m, na.rm=T))/(max(m, na.rm=T) - min(m, na.rm=T))

white2red <- colorRampPalette(c("white", "red"))
cc = c(white2red(500))

image(new.m, xaxt="n", yaxt="n", col = cc)
geneCount = ncol(m)
path.col = rep("black", ncol(m))
mtext(text=colnames(m), side=2, line=0.3, at=(0:(geneCount-1))/(geneCount-1), las=1, cex=.8, col=path.col)
geneCount = nrow(m)

text((0:(geneCount-1))/(geneCount-1), 1.05, srt = 90, adj = 0, cex=.8, labels = rownames(m), xpd = TRUE)

geneCount = nrow(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(v=x, col=grey(0.9), lwd=.2)}
geneCount = ncol(m)
1/(2*geneCount-2) + (0:(geneCount-1))/(geneCount-1) -> tt
for(x in tt){abline(h=x, col=grey(0.8), lwd=.2);}
box()

###############################################################################

m = -log10(TCGA.resistant.mat)

which(is.infinite(m))

m[which( abs(m) < -log10(TCGA.cutoff) )] = 0
m[which( abs(m) > 100 )] = 100
m[which(is.na(m))] = 0
m = t(m)

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

dev.off()

###############################################################################

rbind(cbind(Cancer = names(sensitive.prop), Prop = sensitive.prop, grp="Sensitive"), cbind(Cancer = names(resistant.prop), Prop = -resistant.prop, grp="Insensitive")) -> dat
dat = as.data.frame(dat)
dat[,2] = as.numeric(as.character(dat[,2]))
the_order = names(sensitive.prop)
the_order = the_order[length(the_order):1]
dat[,3] = factor(dat[,3], levels=c("Sensitive", "Insensitive"))

p <- ggplot(dat, aes(x = Cancer, y = Prop, group = grp, fill = grp)) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_x_discrete(limits = the_order) +
  scale_y_continuous(breaks = seq(-0.5, 1, 0.1), labels = abs(seq(-0.5, 1, 0.1))) +
  xlab("") + ylab("") + ggtitle("CCLE model") + 
  theme(legend.position = c(0.9,0.9), legend.title = element_blank(), legend.text=element_text(size=7), 
        plot.title = element_text(hjust = 0.5, size=7),
        panel.background = element_rect(fill =  "grey90"),  axis.text.y=element_text(size=6) ) 

pdf("3D.CCLE.PLX4720.1.pdf", width=4.5, height=3)
print(p)
dev.off()

###############################################################################

pdf("3B.CCLE.PLX4720.2.pdf", width=5, height=3)
h <- hist( PLX4720.pred.dr, breaks=100, plot=FALSE)

cuts <- cut(h$breaks, c(-Inf,quantile(PLX4720.pred.dr, probs=.05),quantile(PLX4720.pred.dr, probs=.95),Inf))
cc = rep( c("skyblue", "grey99", "red"), table(cuts) )
plot(h, col=cc, xlab="PLX4720 predicted ActArea", ylab="", main="Distribution of predicted ActArea\nPLX4720 (CCLE model)", cex.main=.7, cex.axis=.7, cex.lab=.7)

dev.off()

###############################################################################


