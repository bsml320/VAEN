setwd("/work/Figures/Figure6/")

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

minmax_normalization = function(x){(x-min(x))/(max(x)-min(x))}
#################################################################################################

cut.point = 200
dat4boxplot.dat = c()

pdf("KRAS-AZD6244.pdf")
par(mfrow=c(3,1))			

pos.list = list()
pos.list[[ 1 ]] = c(12, 13)
pos.list[[ 2 ]] = c(61)

for(gene in c("KRAS", "HRAS", "NRAS")){
	for(cancer in cancer.types){
		drug = "AZD6244"
		
		pred.dr = cancer.predicted.dr[which(cancer.predicted.dr[,2] == cancer), ]
		
		mutation = read.delim(paste("/work/result.EN/dr.CCLE/04-mix/05.cluster.dSNP/",cancer,".dSNP.4cluster.txt", sep=""), as.is=T)
		gene.mutation = mutation[mutation[,1] == gene, ]
		isoform.mutation = gene.mutation[gene.mutation[,4] == unique(gene.mutation[,4])[1], ]
		isoform.mutation[,7] = unlist(sapply(isoform.mutation[,7], function(u)gsub("p\\.","",u)))
		
		for(posk in c(1,2) ){
			pos = pos.list[[posk]]
			
			ii = which(isoform.mutation$mut_pos >= min(pos) & isoform.mutation$mut_pos <= max(pos))
			if(length(ii) < 5){
				#print(  c(cancer, gene, pos) )
				next
			}
			isoform.mutation = isoform.mutation[ii, ]
			colfunc3 <- colorRampPalette(c("blue", "white", "red"))
			colfunc1 <- colorRampPalette(c("blue", "white"))
			colfunc2 <- colorRampPalette(c("white", "red"))
			
			Y = pred.dr[, drug]
			names(Y) = sapply(pred.dr[,1], function(u)gsub("\\.", "-", u))
			Y = sort(Y)
			Y1 = minmax_normalization(Y) * 100
			cc = c( colfunc1(sum(Y <= mean(Y))), colfunc2(sum(Y > mean(Y))) )
			
			plot(x=1:length(Y), y=rep(0.05, length(Y)), ylim=c(-0.15,0.1), col=cc, type="h", xlim=c(0,length(Y)), yaxt="n", ylab="", xlab="Sample", main=paste(gene, "-", cancer, ", cluster = ", paste(pos, collapse=":"), sep=""))
			match(unique(isoform.mutation[,2]), names(Y)) -> ii
			points(ii, rep(-0.01, length(ii)), col="red", pch=19)
			
			which.min( abs(Y - mean(Y)) ) -> ii
			segments(ii, 0, ii, 0.05, col="purple")
				
			ii = (length(Y) + 1)/2
			segments(ii, 0, ii, 0.05, col="green")
			scaled = length(Y)/(length(Y) - cut.point)
			match(unique(isoform.mutation[,2]), names(Y)) -> ii
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
			
			WT.ss = unique(mutation[-which(mutation[,1] %in% c("KRAS", "NRAS", "BRAF", "HRAS")  ), 2])
			match(  WT.ss, names(Y)) -> WT.ii
			
			if(length(ii) >= 5){
				dat4boxplot.dat = rbind(dat4boxplot.dat, cbind(cancer=cancer, gene=gene, group=paste(pos, collapse=":"), pred.dr = Y[ii]), cbind(cancer=cancer, gene=gene, group="other", pred.dr = Y[WT.ii])  )
			}
		}
	}
}
dev.off()
table(dat4boxplot.dat[,3])



dat4boxplot.dat = as.data.frame(dat4boxplot.dat)
dat4boxplot.dat[,4] = as.numeric(as.character(dat4boxplot.dat[,4]))

apply(dat4boxplot.dat, 1, function(u)paste(u[1], u[2], sep=":")) -> tag
dat4boxplot.dat$tag = tag

p1 = ggplot(dat4boxplot.dat,aes(x=tag ,y=pred.dr,color=group,fill=group))  + 
	  geom_point(position=position_jitterdodge(dodge.width=1), size = 1, stroke = 0, shape = 16) + 
	  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9), lwd=0.4, fatten=0.8) + 
	  xlab("") + ylab("Predicted ActArea") + ggtitle("Ras genes, AZD6244") + 
	  theme(plot.title = element_text(hjust = 0.5, size=8), legend.position = c(0.9, 0.9), axis.text.x = element_text(angle = 30, hjust = 1, size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6) )
      
pdf("6F-AZD6244.pdf", width=5, height=2.5)
print(p1)	  
dev.off()


#################################################################################################

cut.point = 200
dat4boxplot2.dat = c()

pdf("Ras-PD0325901.pdf")
par(mfrow=c(3,1))			

pos.list = list()
pos.list[[ 1 ]] = c(12, 13)
pos.list[[ 2 ]] = c(61)

for(gene in c("KRAS", "HRAS", "NRAS")){
	for(cancer in cancer.types){
		drug = "PD.0325901"
		
		pred.dr = cancer.predicted.dr[which(cancer.predicted.dr[,2] == cancer), ]
		
		mutation = read.delim(paste("/work/result.EN/dr.CCLE/04-mix/05.cluster.dSNP/",cancer,".dSNP.4cluster.txt", sep=""), as.is=T)
		gene.mutation = mutation[mutation[,1] == gene, ]
		isoform.mutation = gene.mutation[gene.mutation[,4] == unique(gene.mutation[,4])[1], ]
		isoform.mutation[,7] = unlist(sapply(isoform.mutation[,7], function(u)gsub("p\\.","",u)))
		
		for(posk in c(1,2) ){
			pos = pos.list[[posk]]
			
			ii = which(isoform.mutation$mut_pos >= min(pos) & isoform.mutation$mut_pos <= max(pos))
			if(length(ii) < 5){
				#print(  c(cancer, gene, pos) )
				next
			}
			isoform.mutation = isoform.mutation[ii, ]
			colfunc3 <- colorRampPalette(c("blue", "white", "red"))
			colfunc1 <- colorRampPalette(c("blue", "white"))
			colfunc2 <- colorRampPalette(c("white", "red"))
			
			Y = pred.dr[, drug]
			names(Y) = sapply(pred.dr[,1], function(u)gsub("\\.", "-", u))
			Y = sort(Y)
			Y1 = minmax_normalization(Y) * 100
			cc = c( colfunc1(sum(Y <= mean(Y))), colfunc2(sum(Y > mean(Y))) )
			
			plot(x=1:length(Y), y=rep(0.05, length(Y)), ylim=c(-0.15,0.1), col=cc, type="h", xlim=c(0,length(Y)), yaxt="n", ylab="", xlab="Sample", main=paste(gene, "-", cancer, ", cluster = ", paste(pos, collapse=":"), sep=""))
			match(unique(isoform.mutation[,2]), names(Y)) -> ii
			points(ii, rep(-0.01, length(ii)), col="red", pch=19)
			
			which.min( abs(Y - mean(Y)) ) -> ii
			segments(ii, 0, ii, 0.05, col="purple")
				
			ii = (length(Y) + 1)/2
			segments(ii, 0, ii, 0.05, col="green")
			scaled = length(Y)/(length(Y) - cut.point)
			match(unique(isoform.mutation[,2]), names(Y)) -> ii
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
			
			WT.ss = unique(mutation[-which(mutation[,1] %in% c("KRAS", "NRAS", "BRAF", "HRAS")  ), 2])
			match(  WT.ss, names(Y)) -> WT.ii
			
			if(length(ii) >= 5){
				dat4boxplot2.dat = rbind(dat4boxplot2.dat, cbind(cancer=cancer, gene=gene, group=paste(pos, collapse=":"), pred.dr = Y[ii]), cbind(cancer=cancer, gene=gene, group="other", pred.dr = Y[WT.ii])  )
			}
		}
	}
}
dev.off()
table(dat4boxplot2.dat[,3])



dat4boxplot2.dat = as.data.frame(dat4boxplot2.dat)
dat4boxplot2.dat[,4] = as.numeric(as.character(dat4boxplot2.dat[,4]))

apply(dat4boxplot2.dat, 1, function(u)paste(u[1], u[2], sep=":")) -> tag
dat4boxplot2.dat$tag = tag

p1 = ggplot(dat4boxplot2.dat,aes(x=tag ,y=pred.dr,color=group,fill=group))  + 
	  geom_point(position=position_jitterdodge(dodge.width=1), size = 1, stroke = 0, shape = 16) + 
	  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9), lwd=0.4, fatten=0.8) + 
	  xlab("") + ylab("Predicted ActArea") + ggtitle("Ras genes, PD0325901") + 
	  theme(plot.title = element_text(hjust = 0.5, size=8), legend.position = c(0.9, 0.9), axis.text.x = element_text(angle = 30, hjust = 1, size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6) )
      
pdf("6F-PD0325901.pdf", width=5, height=2.5)
print(p1)	  
dev.off()

