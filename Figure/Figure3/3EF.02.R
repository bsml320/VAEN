#setwd("/path/to/VAEN/Figure/Figure3")

load("3EF.data.RData")
sapply(rownames(TCGA.sensitive.mat), function(u)gsub("\\.", "-", u)) -> new.name
rownames(TCGA.sensitive.mat) = unlist(new.name)
sapply(rownames(TCGA.resistant.mat), function(u)gsub("\\.", "-", u)) -> new.name
rownames(TCGA.resistant.mat) = unlist(new.name)

two.drugs.match = read.table("../../DATA/drugs.match-2.txt", as.is=T, sep="\t")

gdsc.anno = read.delim("../../DATA/GDSC/Screened_Compounds.txt", as.is=T)
which(gdsc.anno$TARGET_PATHWAY == "ERK MAPK signaling" ) -> ii.1
which(gdsc.anno$TARGET_PATHWAY == "EGFR signaling") -> ii.2
MAPK.inhibitor = sort(unique(gdsc.anno[ii.1, 2]))
EGFR.inhibitor = sort(unique(gdsc.anno[ii.2, 2]))

mek.inhibitor = c(MAPK.inhibitor, EGFR.inhibitor)

mek.inhibitor = c(MAPK.inhibitor, EGFR.inhibitor, setdiff(two.drugs.match[,3], c(MAPK.inhibitor, EGFR.inhibitor) ))
mek.inhibitor[grep("Nutlin", mek.inhibitor)] = rownames(TCGA.sensitive.mat)[grep("Nutlin", rownames(TCGA.sensitive.mat))]

TCGA.cutoff = 0.05/(33*251)
TCGA.cutoff = 0.05

### way 2
match(mek.inhibitor, rownames(TCGA.sensitive.mat)) -> ii
### no overlap
new.TCGA.sensitive.mat = TCGA.sensitive.mat[ii,]
new.TCGA.resistant.mat = TCGA.resistant.mat[ii,]

new.TCGA.sensitive.mat = new.TCGA.sensitive.mat[nrow(new.TCGA.sensitive.mat):1, ]
new.TCGA.resistant.mat = new.TCGA.resistant.mat[nrow(new.TCGA.resistant.mat):1, ]

pdf("3F.GDSC.MEKonly.lineage.pdf", width=5.7, height=5.7)
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


geneCount = ncol(m)
rect(1, 0, 1.01, length(ii.1)/(geneCount-1) , col="red", border=NA)

rect(0.9, 0, 0.91, length(ii.2)/(geneCount-1) , col="lightblue", border=NA)

rect(0.8, 0, 0.81, length( c(ii.2, ii.1) )/(geneCount-1) , col="lightgreen", border=NA)

dev.off()

##############################################################################
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
  xlab("") + ylab("") + ggtitle("GDSC model") + 
  theme(legend.position = c(0.9, 0.9), legend.title = element_blank(), legend.text=element_text(size=7), 
        plot.title = element_text(hjust = 0.5, size=7),
        panel.background = element_rect(fill =  "grey90"),  axis.text.y=element_text(size=6) ) 

pdf("3E.GDSC.PLX4720.pdf", width=4.5, height=3)
print(p)
dev.off()

