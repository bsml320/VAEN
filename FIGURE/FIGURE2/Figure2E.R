setwd("/path/to/VAEN/result.EN/dr.CCLE/VAEN/FIGURE/FIGURE2")
library("Hmisc")

drugs.match = read.table("drugs.match.txt", as.is=T, sep="\t")

ccle = read.table("/path/to/VAEN/result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
gdsc = read.table("/path/to/VAEN/result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")

cancer.types = sort(unique(ccle[,2]))

shared.drugs.mat = matrix(0, nrow=nrow(drugs.match), ncol=length(cancer.types) )
pdf("Figure2E.per.cancer.pdf", width=21, height=15)
for(k1 in 1:nrow(drugs.match)){
	pp.list = list()
	for(k in 1:length(cancer.types)){
		which(ccle[,2]==cancer.types[k]) -> ii;
		
		df = as.data.frame(cbind(x=ccle[ii,  drugs.match[k1,2] ], y=gdsc[ii, drugs.match[k1,3]]) )
		df[,1] = as.numeric(as.character(df[,1]))
		df[,2] = as.numeric(as.character(df[,2]))
		
		rcorr(as.matrix(df), type="pearson") -> a
		shared.drugs.mat[k1, k] = ifelse(a$r[1,2] > 0,  a$P[1,2], 1)
		
		p1 = ggplot(df, aes(x=x, y=y)) + geom_point() + ggtitle(cancer.types[k]) + xlab("CCLE") + ylab("GDSC") +
		   theme(legend.position = "none", text = element_text(size=10), axis.text.x=element_text(color="black", size = 10),axis.text.y=element_text(color="black", size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
		   geom_smooth(method="lm",color='red',data = df, aes(x=x, y=y))
	
		pp.list[[k]] = p1
	}
	
	multiplot(plotlist=pp.list, layout=matrix(c(1:35), ncol=7, byrow=T))
	cat(k1, ".", sep="")
}
dev.off()


shared.drugs.pcc.mat = matrix(0, nrow=nrow(drugs.match), ncol=length(cancer.types) )
for(k1 in 1:nrow(drugs.match)){
	for(k in 1:length(cancer.types)){
		which(ccle[,2]==cancer.types[k]) -> ii;
		
		df = as.data.frame(cbind(x=ccle[ii,  drugs.match[k1,2] ], y=gdsc[ii, drugs.match[k1,3]]) )
		df[,1] = as.numeric(as.character(df[,1]))
		df[,2] = as.numeric(as.character(df[,2]))
		
		rcorr(as.matrix(df), type="pearson") -> a
		shared.drugs.pcc.mat[k1, k] = a$r[1,2] 
	}
	cat(k1, ".", sep="")
}


rownames(shared.drugs.mat) = drugs.match[, 1]
colnames(shared.drugs.mat) = cancer.types
log.shared.drugs.mat = t(-log(shared.drugs.mat + 1e-16))
new = log.shared.drugs.mat[, order(apply(log.shared.drugs.mat, 2, mean))]


library(reshape2)
melt(new) -> dat

pdf("Figure2E.pdf", width=4, height=4) 
ggplot(dat, aes(x=Var2, y=value)) + 
		     geom_boxplot(outlier.shape = NA) + guides(fill=FALSE) +
		     geom_jitter(shape=21, position = position_jitter(width = 0.4), size=1) + guides(colour=FALSE) +
		     xlab("") + ylab("-log(p)") +
		     theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8, angle=90, hjust=1),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(2,2,2,2), "mm"), plot.title = element_text(size = 8, hjust = 0.5))
dev.off()

