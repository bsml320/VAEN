setwd("/work/result.EN/dr.CCLE/04-mix/10.CNV/")
gene.info = read.delim("/work/data/Homo_sapiens.gene_info", as.is=T, skip=1, header=F) ### downloaded from NCBI

##########################################################################################

files = dir("loss.full")
files = files[grep("full.ts.txt", files)]

all.loss.slim.mat = c()
for(kf in 1:length(files)){
	cancer = strsplit(files[kf], split="\\.")[[1]][1]

	loss.sig.mat = read.table(paste("loss.full/",files[kf],sep=""), as.is=T, sep="\t")
	match(loss.sig.mat[,2], gene.info[,3]) -> idx
	loss.sig.mat = cbind(loss.sig.mat, locus = gene.info[idx, 8], type = gene.info[idx, 10])
	loss.sig.mat = loss.sig.mat[which(loss.sig.mat[,19] == "protein-coding"),]

	loss.ps = loss.sig.mat[,4]
	label = paste(loss.sig.mat[,1], loss.sig.mat[,2], loss.sig.mat[,3], sep="-")
	all.loss.slim.mat = rbind(all.loss.slim.mat, cbind(loss.sig.mat, loss.ps, SNP=label) )
}

grep("Un_", all.loss.slim.mat[,10]) -> ii
if(length(ii) > 0)all.loss.slim.mat = all.loss.slim.mat[-ii, ]

gsub("chr", "", all.loss.slim.mat[,10]) -> chr
chr[which(chr=="X")] = 23
all.loss.slim.mat$CHR = as.numeric(chr)


all.loss.slim.mat = all.loss.slim.mat[order(all.loss.slim.mat$CHR, all.loss.slim.mat$V13), ]
d=data.frame(CHR=all.loss.slim.mat[["CHR"]], BP=all.loss.slim.mat[["V13"]], P=all.loss.slim.mat[["loss.ps"]])
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


which(all.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
sort(unique(apply(all.loss.slim.mat[ii, ], 1, function(u)paste(u[1], u[3], sep="-")))) -> tags

library(RColorBrewer)
n <- length(tags)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



tiff("manhantan.loss.sensitive.full.tiff", width=3500, height=1800, res=200)

plot(NA, xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d$logp))), xlab=xlabel, ylab=expression(-log[10](italic(p))))

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
	which(all.loss.slim.mat[,1]==cancer & all.loss.slim.mat[,3] == drug & all.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
	if(length(ii) > 0){
		with(d[ii, ], points(pos, logp, col=col_vector[k], pch=20, cex=.8))
	}
}

icol=1
for (i in unique(d$index)) {
	which( d$index==unique(d$index)[i] )  -> ii
	which.max(d[ii, "logp"]) -> imin
	if(max(d[ii, "logp"]) > 8){
		text.col = "blue"
		if(all.loss.slim.mat[ii[imin], 5] > 0)text.col="red"
		text(d[ii[imin], "pos"], d[ii[imin], "logp"], paste(all.loss.slim.mat[ii[imin], 1:3], collapse="-"), pos=3, cex=.8, col=text.col )
		points(d[ii[imin], "pos"], d[ii[imin], "logp"], cex=.8)
	}
	
	icol=icol+1
}

dev.off()

#########################################################################################

which(all.loss.slim.mat$CHR == 9) -> idx
d3.loss.slim.mat = all.loss.slim.mat[idx, ]
d3=data.frame(CHR=all.loss.slim.mat[idx, "CHR"], BP=all.loss.slim.mat[idx, "V13"], P=all.loss.slim.mat[idx, "loss.ps"])
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

apply(d3.loss.slim.mat, 1, function(u)paste(u[1], u[3], sep="-")) -> d3.tags
d3.tags = unique(d3.tags)

#tiff("chr9.loss.tiff", width=2000, height=1000, res=200)
pdf("chr9.loss.pdf", width=10, height=5)
plot(NA, xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d3$logp))), xlab="Chromosome 9", ylab=expression(-log[10](italic(p))))

col=rep(c("grey", "skyblue"), max(d$CHR))[3]

points(d3$pos, d3$logp, col=col, pch=20)
axis(1, at = c(1000000, 50000000, 100000000, 150000000, 200000000), labels=c("1M", "50M", "100M", "150M", "200M"))

which(d3.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
sort(unique(apply(d3.loss.slim.mat[ii, ], 1, function(u)paste(u[1], u[3], sep="-")))) -> d3.sig.tags

library(RColorBrewer)
n <- length(d3.sig.tags)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


for(k in 1:length(d3.sig.tags)){
	cancer = strsplit(d3.sig.tags[k], split="-")[[1]][1]
	drug   = strsplit(d3.sig.tags[k], split="-")[[1]][2]
	which(d3.loss.slim.mat[,1]==cancer & d3.loss.slim.mat[,3] == drug & d3.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
	
	with(d3[ii, ], points(pos, logp, col=col_vector[k], pch=20))
}

for(k in 1:length(d3.sig.tags)){
	cancer = strsplit(d3.sig.tags[k], split="-")[[1]][1]
	drug   = strsplit(d3.sig.tags[k], split="-")[[1]][2]
	which(d3.loss.slim.mat[,1]==cancer & d3.loss.slim.mat[,3] == drug & d3.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
	which.min(d3.loss.slim.mat[ii, "loss.ps"]) -> imin
	pos1 = 4
	text.col = "blue"
	if(d3.loss.slim.mat[ii[imin], 5] > 0)text.col="red"
	with(d3[ii[imin], ], text(pos, logp, paste(d3.loss.slim.mat[ii[imin], 1:3], collapse="-"), pos=pos1, cex=.6, col=text.col ) )
	with(d3[ii[imin], ], points(pos, logp, col="black", cex=1))
	with(d3[ii[imin], ], points(pos, logp, col=col_vector[k], pch=20))
}

legend(1000000, 25, fill=col_vector, legend=d3.sig.tags, cex=.6)

dev.off()

#########################################################################################

which(all.loss.slim.mat$CHR == 16) -> idx
d3.loss.slim.mat = all.loss.slim.mat[idx, ]
d3=data.frame(CHR=all.loss.slim.mat[idx, "CHR"], BP=all.loss.slim.mat[idx, "V13"], P=all.loss.slim.mat[idx, "loss.ps"])
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

apply(d3.loss.slim.mat, 1, function(u)paste(u[1], u[3], sep="-")) -> d3.tags
d3.tags = unique(d3.tags)

#tiff("chr16.loss.tiff", width=2000, height=1000, res=200)
pdf("chr16.loss.pdf", width=10, height=5)
plot(NA, xaxt='n', bty='n', xaxs='i', yaxs='i', las=1, pch=20, xlim=c(xmin,xmax), ylim=c(0,ceiling(max(d3$logp))), xlab="Chromosome 16", ylab=expression(-log[10](italic(p))))

col=rep(c("grey", "skyblue"), max(d$CHR))[3]

points(d3$pos, d3$logp, col=col, pch=20)
axis(1, at = c(1000000, 25000000, 50000000, 75000000, 100000000), labels=c("1M", "25M", "50M", "75M", "100M"))

which(d3.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
sort(unique(apply(d3.loss.slim.mat[ii, ], 1, function(u)paste(u[1], u[3], sep="-")))) -> d3.sig.tags

library(RColorBrewer)
n <- length(d3.sig.tags)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


for(k in 1:length(d3.sig.tags)){
	cancer = strsplit(d3.sig.tags[k], split="-")[[1]][1]
	drug   = strsplit(d3.sig.tags[k], split="-")[[1]][2]
	which(d3.loss.slim.mat[,1]==cancer & d3.loss.slim.mat[,3] == drug & d3.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
	
	with(d3[ii, ], points(pos, logp, col=col_vector[k], pch=20))
}

for(k in 1:length(d3.sig.tags)){
	cancer = strsplit(d3.sig.tags[k], split="-")[[1]][1]
	drug   = strsplit(d3.sig.tags[k], split="-")[[1]][2]
	which(d3.loss.slim.mat[,1]==cancer & d3.loss.slim.mat[,3] == drug & d3.loss.slim.mat[, "loss.ps"] < 1e-8) -> ii
	which.min(d3.loss.slim.mat[ii, "loss.ps"]) -> imin
	pos1 = 4
	if(d3.sig.tags[k] == "BRCA-PLX4720")pos1 = 3
	text.col = "blue"
	if(d3.loss.slim.mat[ii[imin], 5] > 0)text.col="red"
	with(d3[ii[imin], ], text(pos, logp, paste(d3.loss.slim.mat[ii[imin], 1:3], collapse="-"), pos=pos1, cex=.7, col=text.col ) )
	with(d3[ii[imin], ], points(pos, logp, col="black"))
	with(d3[ii[imin], ], points(pos, logp, col=col_vector[k], pch=20))
}
	
legend(1000000, 15, fill=col_vector, legend=d3.sig.tags, cex=.7)

dev.off()
