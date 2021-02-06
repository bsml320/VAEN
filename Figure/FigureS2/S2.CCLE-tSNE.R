setwd("/path/to/VAEN/Figure/FigureS2")

library("tsne")
library("RColorBrewer")

###################################################################################################

latent = read.table("../../other/RANK.Sigmoid/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)


pdf("RANK.Sigmoid.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white", main="RANK, Sigmoid")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}
plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

###################################################################################################

latent = read.table("../../other/RANK.ReLU/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)


pdf("RANK.ReLU.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white", main="RANK, ReLU")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}
plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

###################################################################################################

latent = read.table("../../other/ZS.Sigmoid/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)


pdf("ZS.Sigmoid.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white", main="ZS, Sigmoid")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}
plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

###################################################################################################

latent = read.table("../../other/ZS.ReLU/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)


pdf("ZS.ReLU.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white", main="ZS, ReLU")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}
plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

###################################################################################################

latent = read.table("../../other/Z01.Sigmoid/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)


pdf("Z01.Sigmoid.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white", main="Z01, Sigmoid")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}
plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

###################################################################################################

latent = read.table("../../other/Z01.ReLU/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)


pdf("Z01.ReLU.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white", main="Z01, ReLU")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}
plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

###################################################################################################
