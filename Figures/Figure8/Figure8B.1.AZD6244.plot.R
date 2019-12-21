setwd("/work/result.EN/dr.CCLE/04-mix/09.Expr/AZD6244")

mydata = read.table("/work/result.EN/dr.CCLE/04-mix/09.Expr/W5/AZD6244.txt", header=T, as.is=T)
colnames(mydata)[1:3] = c("cancer", "gene", "drug")
mydata = mydata[mydata[,1]!="BRCA" & mydata[,1]!="LGG" & mydata[,1]!="THCA",]

tapply(mydata[,5], mydata[,"tag"], function(u)sum(u)/sqrt(length(u))) -> check
check = sort(check)
print(tail(check))

AZD6244.18 = read.table("/work/result.EN/dr.CCLE/04-mix/09.Expr/AZD6244/AZD6244-18.PMC3166660.txt", as.is=T)
AZD6244.18 = AZD6244.18[,1]

AZD6244.sen = read.table("/work/result.EN/dr.CCLE/04-mix/09.Expr/AZD6244/AZD6244-sensitive.PMC3931013.txt", as.is=T)
AZD6244.res = read.table("/work/result.EN/dr.CCLE/04-mix/09.Expr/AZD6244/AZD6244-resistent.PMC3931013.txt", as.is=T)
AZD6244.sen = AZD6244.sen[,1]
AZD6244.res = AZD6244.res[,1]

pdf("Figure8B.1.AZD6244.pdf", width=5, height=5)
plot(check, type="l", xlab="Gene index", ylab="Average t-value", main="AZD6244")
match(paste(AZD6244.18, "_AZD6244", sep=""), names(check)) -> ii
ii = ii[!is.na(ii)]

for(k in ii){
	segments(k, 20, k, 21, col="green")
}

ii/length(check)

match(paste(AZD6244.sen, "_AZD6244", sep=""), names(check)) -> ii
ii = ii[!is.na(ii)]
for(k in ii){
	segments(k, 18, k, 19, col="red")
}

ii/length(check)

match(paste(AZD6244.res, "_AZD6244", sep=""), names(check)) -> ii
ii = ii[!is.na(ii)]

for(k in ii){
	segments(k, 16, k, 17, col="purple")
}

ii/length(check)

match("SPRY2_AZD6244", names(check)) -> k
points(k, check[k], pch=19, col="red")
text(k, check[k], "SPRY2", cex=.8)

dev.off()

#######################################################################################
