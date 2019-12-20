setwd("/work/Figures/Figure6/")
minmax_normalization = function(x, new.min, new.max){(x-min(x))/(max(x)-min(x)) -> base; base * (new.max-new.min) + new.min}

predicted.dr = read.table("/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
drugs = colnames(predicted.dr)[c(-1, -2)]
cancer.types = unique(predicted.dr[,2])
sample.type = substr(predicted.dr[,1], 14, 15)

cancer.predicted.dr = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = predicted.dr[which(predicted.dr[,2] == cancer & sample.type == type.code), ]
	cancer.predicted.dr = rbind(cancer.predicted.dr, blca.ccle)
}

############################################################################

colfunc1 <- colorRampPalette(c("blue", "white"))
colfunc2 <- colorRampPalette(c("white","red"))
colfunc3 <- colorRampPalette(c("blue", "white", "red"))

BLCA.dr = cancer.predicted.dr[which(cancer.predicted.dr[,2] == "BLCA"), ]
res.mat = read.table("/work/result.EN/dr.CCLE/04-mix/05.cluster.dSNP/04.05.cluster.dSNP.txt", as.is=T)
FGFR3.sig.mat = res.mat[which(res.mat[,3]=="NM_000142"),]

pdf("BLCA-FGFR3.pdf", width=8, height=5)
	
	for(krow in 1:nrow(FGFR3.sig.mat)){
		Y = BLCA.dr[, FGFR3.sig.mat[krow,4]]
		names(Y) = sapply(BLCA.dr[,1], function(u)gsub("\\.", "-", u))
		Y = sort(Y)
		Y1 = minmax_normalization(Y, 0, 1) * 100
		cc = colfunc3(101)[ 1+ceiling(Y1) ]
		
		cc = c( colfunc1(sum(Y <= mean(Y))), colfunc2(sum(Y > mean(Y))) )
		
		
		mutation = read.delim("/work/result.EN/dr.CCLE/04-mix/05.cluster.dSNP/BLCA.dSNP.4cluster.txt", as.is=T)  ### output from 04.05.cluster.dSNP.R
		isoform.mutation = mutation[mutation[,4] == FGFR3.sig.mat[krow, 3], ]
		isoform.mutation[,7] = unlist(sapply(isoform.mutation[,7], function(u)gsub("p\\.","",u)))
		
		strsplit(FGFR3.sig.mat[krow,9], split=":")[[1]] -> pos
		isoform.mutation = isoform.mutation[isoform.mutation[,8] %in% pos, ]
		
		plot(x=1:length(Y), y=rep(0.05, length(Y)), ylim=c(-0.15,0.1), col=cc, type="h", xlim=c(0,length(Y)), yaxt="n", ylab="", xlab="Sample", main=paste(FGFR3.sig.mat[krow, 2], "-", FGFR3.sig.mat[krow, 4], ", cluster = ", FGFR3.sig.mat[krow, 9], sep=""))
		
		which.min( abs(Y - mean(Y)) ) -> ii
		segments(ii, 0, ii, 0.05, col="purple")
		
		ii = (length(Y) + 1)/2
		segments(ii, 0, ii, 0.05, col="green")
		
		match(unique(isoform.mutation[,2]), names(Y)) -> ii
		points(ii, rep(-0.005, length(ii)), col="red", pch=19)
		
		cut.point = 200
		segments(cut.point,-0.045,cut.point, -0.012)
		segments(1,-0.08,cut.point, -0.045)
		segments(length(Y),-0.08,length(Y), -0.012)
		segments(cut.point, -0.012, length(Y), -0.012)
		
		scaled = length(Y)/(length(Y) - cut.point)
		
		if(length(ii) < 20){
			for(k in ii){
				which(isoform.mutation[, 2] == names(Y)[k]) -> mut.ii
				text( k, -0.03, labels=isoform.mutation[mut.ii[1], 7], srt=90, cex=.8)
			}
		
		} else {
		
		for(k in ii){
			if(k > cut.point){
				which(isoform.mutation[, 2] == names(Y)[k]) -> mut.ii
				text( (k - cut.point) * scaled, -0.1, labels=isoform.mutation[mut.ii[1], 7], srt=90, cex=.8)
				points((k - cut.point) * scaled, -0.08, pch=19, col="red")
			} else {
				which(isoform.mutation[, 2] == names(Y)[k]) -> mut.ii
				text( k, -0.03, labels=isoform.mutation[mut.ii[1], 7], srt=90, cex=.8)
			}
		}
		
		}
	}
dev.off()
