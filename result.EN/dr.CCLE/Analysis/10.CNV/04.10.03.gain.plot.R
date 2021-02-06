setwd("/path/to/VAEN/result.EN/dr.CCLE/Analysis/10.CNV")
source("/path/to/VAEN/code/unfactor.R")

gene.info = read.delim("/path/to/VAEN/DATA/Homo_sapiens.gene_info", as.is=T, skip=1, header=F)

##########################################################################################

files = dir("gain")
files = files[grep("full.ts.txt", files)]

all.gain.slim.mat = c()
for(kf in 1:length(files)){
	cancer = strsplit(files[kf], split="\\.")[[1]][1]

	gain.sig.mat = read.table(paste("gain/",files[kf],sep=""), as.is=T, sep="\t")
	match(gain.sig.mat[,2], gene.info[,3]) -> idx
	gain.sig.mat = cbind(gain.sig.mat, locus = gene.info[idx, 8], type = gene.info[idx, 10])
	gain.sig.mat = gain.sig.mat[which(gain.sig.mat[,19] == "protein-coding"),]

	gain.ps = gain.sig.mat[,4]
	label = paste(gain.sig.mat[,1], gain.sig.mat[,2], gain.sig.mat[,3], sep="-")
	all.gain.slim.mat = rbind(all.gain.slim.mat, cbind(gain.sig.mat, gain.ps, SNP=label) )
}

gsub("chr", "", all.gain.slim.mat[,10]) -> chr
chr[which(chr=="X")] = 23
all.gain.slim.mat$CHR = as.numeric(chr)


all.gain.slim.mat = all.gain.slim.mat[order(all.gain.slim.mat$CHR, all.gain.slim.mat$V13), ]
d=data.frame(CHR=all.gain.slim.mat[["CHR"]], BP=all.gain.slim.mat[["V13"]], P=all.gain.slim.mat[["gain.ps"]])
d$logp <- -log10(d$P)
d$index=NA
ind = 0
for (i in unique(d$CHR)){
	ind = ind + 1
	d[d$CHR==i,]$index = ind
}


nchr = length(unique(d$CHR))
lastbase=0
ticks=NULL
d$pos = 0
for (i in unique(d$index)) {
	if (i==1) {
        d[d$index==i, "pos"]=d[d$index==i, "BP"]
    } else {
        lastbase=lastbase+tail(subset(d,index==i-1)$BP, 1)
        d[d$index==i, "pos"] = d[d$index==i, "BP"]+lastbase
    }
    ticks = c(ticks, (min(d[d$index == i,]$pos) + max(d[d$index == i,]$pos))/2 + 1)
}

xlabel = 'Chromosome'
labs <- unique(d$CHR)

xmax = ceiling(max(d$pos) * 1.03)
xmin = floor(max(d$pos) * -0.03)


which(all.gain.slim.mat[, "gain.ps"] < 1e-8) -> ii
sort(unique(apply(all.gain.slim.mat[ii, ], 1, function(u)paste(u[1], u[3], sep="-")))) -> tags

library(RColorBrewer)
n <- length(tags)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



tiff("7E.tiff", width=2000, height=1000, res=200)
par(mar=c(5.1, 4.1, 1.1, 0.1))
#plot(NA, xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$logp))), xlab=xlabel, ylab=expression(-log[10](italic(p))))
plot(NA, xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(0,30), xlab=xlabel, ylab=expression(-log[10](italic(p))))
col=rep(c("grey", "skyblue"), max(d$CHR))
icol=1
for (i in unique(d$index)) {
	with(d[d$index==unique(d$index)[i] & d$logp < 8, ], points(pos, logp, col=col[icol], pch=20))
	icol=icol+1
}
axis(1, at=ticks, labels=labs)

for(k in 1:length(tags)){
	cancer = strsplit(tags[k], split="-")[[1]][1]
	drug = strsplit(tags[k], split="-")[[1]][2]
	which(all.gain.slim.mat[,1]==cancer & all.gain.slim.mat[,3] == drug & all.gain.slim.mat[, "gain.ps"] < 1e-10) -> ii
	if(length(ii) > 0){
		with(d[ii, ], points(pos, logp, col=col_vector[k], pch=20, cex=.8))
	}
}

icol=1
for (i in unique(d$index)) {
	which( d$index==unique(d$index)[i] )  -> ii
	which.max(d[ii, "logp"]) -> imin
	if(max(d[ii, "logp"]) > 10){
		text.col = "blue"
		if(all.gain.slim.mat[ii[imin], 5] > 0)text.col="red"
		if(i == 17 & all.gain.slim.mat[ii[imin], 21]=="BRCA-IKZF3-Lapatinib")next
		text(d[ii[imin], "pos"], d[ii[imin], "logp"], paste(all.gain.slim.mat[ii[imin], 1:3], collapse="-"), pos=3, cex=.8, col=text.col )
		points(d[ii[imin], "pos"], d[ii[imin], "logp"], cex=.8)
	}
	
	icol=icol+1
}

which(all.gain.slim.mat$SNP=="BRCA-ERBB2-Lapatinib") -> ii
text(d[ii, "pos"], d[ii, "logp"], paste(all.gain.slim.mat[ii, 1:3], collapse="-"), pos=3, cex=.8, col="red" )
points(d[ii, "pos"], d[ii[imin], "logp"], cex=.8)

dev.off()

#########################################################################################
#########################################################################################

which(all.gain.slim.mat$CHR == 17) -> idx
d3.gain.slim.mat = all.gain.slim.mat[idx, ]
d3=data.frame(CHR=all.gain.slim.mat[idx, "CHR"], BP=all.gain.slim.mat[idx, "V13"], P=all.gain.slim.mat[idx, "gain.ps"])
d3$logp <- -log10(d3$P)
d3$index=NA
d3$pos = d3$BP
ind = 0
for (i in unique(d3$CHR)){
	ind = ind + 1
	d3[d3$CHR==i,]$index = ind
}

xmin = 1
xmax = max(d3$pos) * 1.05

apply(d3.gain.slim.mat, 1, function(u)paste(u[1], u[3], sep="-")) -> d3.tags
d3.tags = unique(d3.tags)

tiff("7F.chr17.tiff", width=2000, height=1000, res=200)
par(mar=c(5.1, 4.1, 1.1, 0.1))
plot(NA, xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d3$logp))), xlab="Chromosome 17", ylab=expression(-log[10](italic(p))))

col=rep(c("grey", "skyblue"), max(d$CHR))[3]

points(d3$pos, d3$logp, col=col, pch=20)
axis(1, at = c(1000000, 25000000, 50000000, 75000000, 100000000), labels=c("1M", "25M", "50M", "75M", "100M"))

which(d3.gain.slim.mat[, "gain.ps"] < 1e-8) -> ii
sort(unique(apply(d3.gain.slim.mat[ii, ], 1, function(u)paste(u[1], u[3], sep="-")))) -> d3.sig.tags

library(RColorBrewer)
n <- length(d3.sig.tags)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


for(k in 1:length(d3.sig.tags)){
	cancer = strsplit(d3.sig.tags[k], split="-")[[1]][1]
	drug   = strsplit(d3.sig.tags[k], split="-")[[1]][2]
	which(d3.gain.slim.mat[,1]==cancer & d3.gain.slim.mat[,3] == drug & d3.gain.slim.mat[, "gain.ps"] < 1e-8) -> ii
	
	with(d3[ii, ], points(pos, logp, col=col_vector[k], pch=20))
}

for(k in 1:length(d3.sig.tags)){
	cancer = strsplit(d3.sig.tags[k], split="-")[[1]][1]
	drug   = strsplit(d3.sig.tags[k], split="-")[[1]][2]
	which(d3.gain.slim.mat[,1]==cancer & d3.gain.slim.mat[,3] == drug & d3.gain.slim.mat[, "gain.ps"] < 1e-8) -> ii
	which.min(d3.gain.slim.mat[ii, "gain.ps"]) -> imin
	pos1 = 4
	if(d3.sig.tags[k] == "BRCA-PLX4720")pos1 = 3
	text.col = "blue"
	if(d3.gain.slim.mat[ii[imin], 5] > 0)text.col="red"
	with(d3[ii[imin], ], text(pos, logp, paste(d3.gain.slim.mat[ii[imin], 1:3], collapse="-"), pos=pos1, cex=.7, col=text.col ) )
	with(d3[ii[imin], ], points(pos, logp, col="black"))
	with(d3[ii[imin], ], points(pos, logp, col=col_vector[k], pch=20))
}

which(d3.gain.slim.mat[,1]=="BRCA" & d3.gain.slim.mat[,3] == "Lapatinib" & d3.gain.slim.mat[, 2] == "ERBB2") -> ii
text.col = "blue"
if(d3.gain.slim.mat[ii, 5] > 0)text.col="red"
with(d3[ii, ], text(pos, logp, paste(d3.gain.slim.mat[ii, 1:3], collapse="-"), pos=1, cex=.7, col=text.col ) )
with(d3[ii, ], points(pos, logp, col="red", pch=17))
	
legend(1000000, 15, fill=col_vector, legend=d3.sig.tags, cex=.7)

dev.off()
