setwd("/path/to/VAEN/Figure/Figure2")
source("../../code/unfactor.R")

######################### 2A
one.drugs.match = read.table("../../DATA/drugs.match.txt", as.is=T)
two.drugs.match = read.table("../../DATA/drugs.match-2.txt", as.is=T, sep="\t")

load("../../result.EN/dr.CCLE/dr.CCLE.A.models.RData")
drugs = names(dr.ccle.models)
x = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	res.list = dr.ccle.models[[ drug ]]
	Ys = res.list$Ys
	which(Ys[,1]!=-9) -> ii
	Ys = Ys[ii, ]
	recall    = cor(Ys[,1], Ys[,2])
	
	x = rbind(x, c(drug, recall) )
}
colnames(x) = c("Drug", "PCC")


load("../../result.EN/dr.GDSC/dr.GDSC.A.models.RData")
drugs = names(dr.gdsc.models)
y = c()
for(k in 1:length(drugs)){
	drug = drugs[k]
	res.list = dr.gdsc.models[[ drug ]]
	Ys = res.list$Ys
	which(Ys[,1]!=-9) -> ii
	Ys = Ys[ii, ]
	recall    = cor(Ys[,1], Ys[,2])
	
	y = rbind(y, c(drug, recall) )
}
colnames(y) = c("Drug", "PCC")

write.table(x, file="Figure2A.CCLE.txt")
write.table(y, file="Figure2A.GDSC.txt")

pdf("2A.PCC.pdf", width=8, height=4)
match(two.drugs.match[,1], x[,1]) -> sx.ii
match(two.drugs.match[,3], y[,1]) -> sy.ii

a1 = length(x[-sx.ii, "PCC"])
b1 = length(y[-sy.ii, "PCC"])
a2 = length(sx.ii)
b2 = length(sy.ii)

dat.plot = rbind(
                 cbind(x[-sx.ii, c("Drug", "PCC")], grp="1"), 
				 cbind(y[-sy.ii, c("Drug", "PCC")], grp="2"), 
				 cbind(x[sx.ii, c("Drug", "PCC")], grp="3"), 
				 cbind(y[sy.ii, c("Drug", "PCC")], grp="4")
				)
dat.plot = unfactor(as.data.frame(dat.plot))
dat.plot = dat.plot[order(dat.plot[,3], dat.plot[,2]),]

plot(0,0, xlim=c(-1, a1+b1+a2+b2+10), ylim=c(0,1.05), col="white", ylab="", xlab="")
mtext("Compound index", 1, line=2)
mtext("In-sample PCC", 2, line=2)

rect(-1, 0, a1+0.5, 1, col="lightyellow", border=F)
rect(a1+0.5, 0, a1+b1+0.5, 1, col="lightcyan", border=F)
rect(a1+b1+0.5, 0, a1+b1+a2+0.5, 1, col="lightgreen", border=F)
rect(a1+b1+a2+0.5, 0, a1+b1+a2+b2+1, 1, col="lightblue", border=F)

points(dat.plot[,2], pch=19, col="grey", cex=1)
points(dat.plot[,2], col="black", cex=1)
segments(-1,0.5,276,0.5)

which(dat.plot[,2] < 0.5 & dat.plot[,3] == "1") -> ii
for(k in ii)text(k, as.numeric(dat.plot[k,2]), dat.plot[k,1], pos=4, cex=.8)

which(dat.plot[,2] < 0.5 & dat.plot[,3] == "3") -> ii
text(ii, as.numeric(dat.plot[ii,2]), dat.plot[ii,1], pos=4, cex=.8)

which(dat.plot[,2] < 0.5 & dat.plot[,3] == "4") -> ii
text(ii, as.numeric(dat.plot[ii,2]), dat.plot[ii,1], pos=4, cex=.8)

dev.off()
