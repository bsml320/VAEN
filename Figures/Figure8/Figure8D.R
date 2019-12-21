setwd("/work/result.EN/dr.CCLE/04-mix/09.Expr/")
W5.files = dir("/work/result.EN/dr.CCLE/04-mix/09.Expr/W5/")

TIS.genes = read.table("/work/data/TIS18.weight.txt")
TIS.genes = as.character(TIS.genes[,1])

################################################################################## mean rank
mean.check = c()
for(k in 1:length(W5.files)){
	drug = gsub("W5\\.", "", gsub(".txt", "", W5.files[k]))
	mydata = read.table(paste("/work/result.EN/dr.CCLE/04-mix/09.Expr/W5/", W5.files[k], sep=""), header=T, as.is=T)
	colnames(mydata)[1:3] = c("cancer", "gene", "drug")
	
	mydata = mydata[mydata[,1]!="BRCA" & mydata[,1]!="LGG" & mydata[,1]!="THCA",]
	tapply(mydata[,"tvalue"], mydata[,"gene"], function(u)sum(u)/sqrt(length(u))) -> check
	check = sort(check)
	mean.check = rbind(mean.check, cbind(gene=names(check), check = check, drug=drug) )
}

mean.check = as.data.frame(mean.check)
source("C:/Users/pjia/UTH/code/unfactor.R")
mean.check = unfactor(mean.check)

################################################################################################### mean data

pdf("8D.pdf", height=7)
drugs = unique(mean.check[,3])

plot(0, 0, xlim=range(mean.check[,2]), ylim=c(0, length(drugs)+1 ), main="Average TIS t-value, rank", xlab="", ylab="", col="white")
	for(kdrug in 1:length(drugs)){
		drug = drugs[kdrug]
		cancer.mydata = mean.check[which(mean.check[,3]==drug), ]
		cancer.mydata = cancer.mydata[order(as.numeric(cancer.mydata[, 2])), ]
		
		cancer.t = as.numeric(cancer.mydata[, 2])
		names(cancer.t) = cancer.mydata[,1]
		
		segments( -1, kdrug, 1, kdrug, col="grey")
		segments( min(cancer.t), kdrug, 0, kdrug, lwd=3, col="lightgreen")
		segments( 0, kdrug, max(cancer.t), kdrug, lwd=3, col="pink")
		
		which(cancer.mydata[,1] %in% TIS.genes) -> ii
		if(length(ii) < 1)next
		for(point.i in ii){
			point1.x = point.i/length(cancer.t) * (max(cancer.t) - min(cancer.t) ) + min(cancer.t)
			point1.y = kdrug
			segments( point1.x, kdrug-0.2, point1.x, kdrug+0.2, lwd=1, col="blue")
		}
		text( min(mean.check[,2]), kdrug, drug, pos=2 )
		
		if(sum(cancer.mydata[ii, 2] < 0) <= 3){
			for(point.i in ii){
				if(cancer.mydata[point.i, 2] < 0){
					point1.x = point.i/length(cancer.t) * (max(cancer.t) - min(cancer.t) ) + min(cancer.t)
					point1.y = kdrug
					text( point1.x, point1.y, cancer.mydata[point.i, 1], pos=2, cex=.8 )
				}
			}
		}
	}
	abline(v=0, col="lightblue")

dev.off()
