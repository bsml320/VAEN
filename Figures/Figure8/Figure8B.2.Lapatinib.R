setwd("/work/result.EN/dr.CCLE/04-mix/09.Expr/Lapatinib")

mydata = read.table("/work/result.EN/dr.CCLE/04-mix/09.Expr/W5/Lapatinib.txt", header=T, as.is=T)
colnames(mydata)[1:3] = c("cancer", "gene", "drug")
mydata = mydata[mydata[,1]!="BRCA" & mydata[,1]!="LGG" & mydata[,1]!="THCA",]

tapply(mydata[,5], mydata[,"tag"], function(u)sum(u)/sqrt(length(u))) -> check
check = sort(check)
print(head(check))


STROMA.genes = read.table("/work/result.EN/dr.CCLE/04-mix/09.Expr/Lapatinib/STROMA.txt", as.is=T)
STROMA.genes = STROMA.genes[, 1]


pdf("Lapatinib.pdf", width=5, height=5)
plot(check, type="l", xlab="Gene index", ylab="Average t-value", main="Lapatinib, CAF")
match(paste(caf.genes, "_Lapatinib", sep=""), names(check)) -> ii
ii = ii[!is.na(ii)]
points(ii, check[ii], pch=19, col="red")

for(k in ii){
	segments(k, 20, k, 21, col="red")
	text(k, check[k], strsplit(names(check)[k], split="_")[[1]][1], cex=.8, pos=4)
}


plot(check, type="l", xlab="Gene index", ylab="Average t-value", main="Lapatinib, STROMA")
match(paste(STROMA.genes, "_Lapatinib", sep=""), names(check)) -> ii
ii = ii[!is.na(ii)]
points(ii, check[ii], pch=19, col="red")

for(k in ii){
	segments(k, 20, k, 21, col="red")
	text(k, check[k], strsplit(names(check)[k], split="_")[[1]][1], cex=.8, pos=4)
}

dev.off()


##############################################################################################################
