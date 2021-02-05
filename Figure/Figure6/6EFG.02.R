setwd("/path/to/VAEN/Figure/Figure6")

dat.mat = read.table("6E.GDSC.tml.ttest.txt", header=T, as.is=T)
dat.mat = dat.mat[-which(dat.mat[,1] %in% c("BRCA", "THCA", "LGG") ),]

dat.mat$adjp = 1

ss = dat.mat[,2]
gsub("\\.", "-", ss) -> new.ss
dat.mat[,2] = new.ss

####### *****
dat.mat$adjp = NULL
drugs = unique(dat.mat[,1])
for(k in 1:length(drugs)){
	which(dat.mat[,1] == drugs[k]) -> ii
	p.adjust(dat.mat[ii,3], method="BH") -> adjp
	dat.mat[ii, "adjp"] = adjp
}
####### *****

which(dat.mat[,2]=="X5.Fluorouracil") -> ii; dat.mat[ii, 2] = "5-Fluorouracil"
which(dat.mat[,2]=="X.5Z..7.Oxozeaenol") -> ii; dat.mat[ii, 2] = "(5Z)-7-Oxozeaenol"
which(dat.mat[,2]=="JW.7.52.1") -> ii; dat.mat[ii, 2] = "JW-7-52-1"


which(dat.mat[,2]=="X-5Z--7-Oxozeaenol") -> ii; dat.mat[ii, 2] = "(5Z)-7-Oxozeaenol"
which(dat.mat[,2]=="X5-Fluorouracil") -> ii; dat.mat[ii, 2] = "5-Fluorouracil"
which(dat.mat[,2]=="Nutlin-3a----") -> ii; dat.mat[ii, 2] = "Nutlin-3a (-)"
which(dat.mat[,2]=="Bleomycin--50-uM-") -> ii; dat.mat[ii, 2] = "Bleomycin (50 uM)"
which(dat.mat[,2]=="Vinorelbine") -> ii; dat.mat[ii, 2] = "Vinorelbine"


pdf("6E.pdf")
plot(x=log2(dat.mat[,4]), y=-log10(dat.mat[,3]), col="#00AFBB", xlab="", ylab="-log10(p)")
abline(v=0, lty=2, col="red")
mtext("log2(FC)",1,line=2)

which(dat.mat[,3] < 1e-9 & dat.mat[,4] > 1) -> ii
text(log2(dat.mat[ii, 4]), -log10(dat.mat[ii, 3]), labels=dat.mat[ii, 2], pos=4, cex=.9)

which(dat.mat[,3] < 1e-19 & dat.mat[,4] < 1) -> ii
text(log2(dat.mat[ii, 4]), -log10(dat.mat[ii, 3]), labels=dat.mat[ii, 2], pos=4, cex=.9)

which(log2(dat.mat[,4]) > 0.2 & dat.mat[,5] < 0.05)  -> ii
points(x=log2(dat.mat[ii,4]), y=-log10(dat.mat[ii,3]), col="red", pch=19)
which(log2(dat.mat[,4]) < -0.2 & dat.mat[,5] < 0.05)  -> ii
points(x=log2(dat.mat[ii,4]), y=-log10(dat.mat[ii,3]), col="blue", pch=19)

dev.off()

###################################################################################################

gdsc.anno = read.delim("../../DATA/GDSC/Screened_Compounds-2.txt", as.is=T)
match(dat.mat[,2], gdsc.anno[,2]) -> ii
cbind(dat.mat, gdsc.anno[ii,2], gdsc.anno[ii, ]) -> new.mat


pdf("6G.pdf", width=3.6, height=4)
par(mfrow=c(1, 1), mar=c(0.2,3,1,1))
par(lheight=.8)
which(log2(dat.mat[,4]) > 0.2 & dat.mat[,3] < 0.01)  -> ii
print(length(ii))
table(new.mat[ii, 11]) -> a
a = a[a>1]
tag = names(a)
res = c()
for(k in 1:length(tag)){
	sum(new.mat[ii,11]==tag[k], na.rm=T)/length(ii) -> r1
	sum(new.mat[,11]==tag[k], na.rm=T)/nrow(new.mat) -> r2
	
	x = ifelse(log2(dat.mat[,4]) > 0.2 & dat.mat[,3] < 0.01, 0, 1 )
	y = ifelse(new.mat[,11]==tag[k], 0, 1)
	
	if( sum(gdsc.anno[,5] == tag[k]) < 5)next
	
	table(x, y) -> mm
	fisher.test(mm)$p.value -> p
	
	res = rbind(res, c(tag[k], r1, r2, r1/r2, p, as.vector(mm) ))
}
res = res[order(as.numeric(res[,4])),]
left = sum(as.numeric(res[,4]) < 1)
right = sum(as.numeric(res[,4]) > 1)
barplot(as.numeric(res[,4]) - 1, ylim=c(-1.5,1.5), col=c(rep("lightblue",sum(as.numeric(res[,4]) < 1)), rep("red", sum(as.numeric(res[,4]) > 1) ) ), ylab="") -> h
#res[which(res[,1]=="Protein stability and degradation"), 1] = "Protein stability\nand degradation"
#res[which(res[,1]=="Chromatin histone methylation"), 1] = "Chromatin histone\nmethylation"
#res[which(res[,1]=="Chromatin histone acetylation"), 1] = "Chromatin histone\nacetylation"

res[which(res[,1]=="JNK and p38 signaling"), 1] = "JNK & p38 signaling"
text(h[(left+1):nrow(res)], -0.1, label=res[(left+1):nrow(res),1], srt=90, adj=1, cex=.8)
text(h[1:left], 0.1, label=res[1:left,1], srt=90, adj=0, cex=.8)
mtext("Odds Ratio - 1",2,line=2)

which(as.numeric(res[,5]) < 0.05) -> ii
if(length(ii) > 0)for(k in ii)text(h[k, 1],as.numeric(res[k,4]) - 1, paste("p = ",format(as.numeric(res[k,5]), digits=3), sep=""), srt=90 )
mtext("(G)", side=2, at = 1.2, las=1, line=2, cex=.9)

dev.off()

########### plot2
pdf("6F.pdf", width=4, height=4)
par(mfrow=c(1, 1), mar=c(0.2,3,1,1))
par(lheight=.8)

which(log2(dat.mat[,4]) < -0.2 & dat.mat[,3] < 0.01)  -> ii
print(length(ii))
table(new.mat[ii, 11]) -> a
a = a[a>1]
tag = names(a)
res2 = c()
for(k in 1:length(tag)){
	sum(new.mat[ii,11]==tag[k], na.rm=T)/length(ii) -> r1
	sum(new.mat[,11]==tag[k], na.rm=T)/nrow(new.mat) -> r2
	
	x = ifelse(log2(dat.mat[,4]) < -0.2 & dat.mat[,3] < 0.01, 0, 1 )
	y = ifelse(new.mat[,11]==tag[k], 0, 1)
	table(x, y) -> mm
	fisher.test(mm)$p.value -> p
	
	res2 = rbind(res2, c(tag[k], r1, r2, r1/r2, p, as.vector(mm) ))
}
res2 = res2[order(as.numeric(res2[,4])),]
left = sum(as.numeric(res2[,4]) < 1)
right = sum(as.numeric(res2[,4]) > 1)
barplot(as.numeric(res2[,4]) - 1, ylim=c(-1.5,1.5), col=c(rep("lightblue",left), rep("red", right ) ), ylab="") -> h
#res2[which(res2[,1]=="Chromatin histone acetylation"), 1] = "Chromatin histone\nacetylation"
res2[which(res2[,1]=="JNK and p38 signaling"), 1] = "JNK & p38 signaling"
#res2[which(res2[,1]=="Protein stability and degradation"), 1] = "Protein stability\nand degradation"
text(h[(left+1):nrow(res2)], -0.1, label=res2[(left+1):nrow(res2),1], srt=90, adj=1, cex=.8)
text(h[1:left], 0.1, label=res2[1:left,1], srt=90, adj=0, cex=.8)
mtext("Odds Ratio - 1",2,line=2)

which(as.numeric(res2[,5]) < 0.05) -> ii
if(length(ii) > 0)for(k in ii)text(h[k, 1],as.numeric(res2[k,4]) - 1, paste("p = ",format(as.numeric(res2[k,5]), digits=3), sep=""), srt=90 )

mtext("(F)", side=2, at = 1.1, las=1, line=2, cex=.9)
mtext("(A)", side=2, at = 1.0, las=1, line=2, cex=.9)
mtext("(B)", side=2, at = 1.3, las=1, line=2, cex=.9)
mtext("(C)", side=2, at = 1.4, las=1, line=2, cex=.9)
mtext("(D)", side=2, at = 1.5, las=1, line=2, cex=.9)
mtext("(E)", side=2, at = 1.6, las=1, line=2, cex=.9)

dev.off()



