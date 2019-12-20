setwd("/work/Figures/Figure5")
library("RColorBrewer")
library("reshape2")
library("gplots")
library("ggplot2")

### RUN /work/result.EN/dr.CCLE/04-mix/12.TML/04.12.CCLE.TumorMutationLoad.R

load("/work/result.EN/dr.CCLE/04-mix/12.TML/cancer.tml.list.RData")

##############################
drug.ccle = read.table("/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T)
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss


pdf("Figure5CD.ccle.drug.pdf", width=5, height=4)
for(kdrug in 1:length(drugs)){
	dat4plot = c()
	for(ct in 1:length(cancer.types)){
		cancer = cancer.types[ct]
		sample2gene.length = cancer.tml.list[[cancer]]
		
		fixed.ss = names(sample2gene.length)
		if(length(fixed.ss) < 50){
			next
		}
		
		blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
		apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle
		
		sample2gene.length.log = log(sample2gene.length)
		upper.border = mean(sample2gene.length.log) + 3 * sd(sample2gene.length.log)
		which(sample2gene.length.log > upper.border) -> ii
		sample2gene.length.log[ii] = upper.border
			
		which(new.ccle[, kdrug] > quantile(new.ccle[, kdrug], probs=.9)) -> ii
		dat4plot = rbind(dat4plot, cbind(cancer, drug=colnames(new.ccle)[kdrug], group = "responder", value=sample2gene.length.log[ii] ) )
		dat4plot = rbind(dat4plot, cbind(cancer, drug=colnames(new.ccle)[kdrug], group = "nonresponder", value=sample2gene.length.log[-ii]  ) )
	}
	dat4plot = dat4plot[!is.na(dat4plot[,4]), ]
	dat4plot = as.data.frame(dat4plot)
	dat4plot[,4] = as.numeric(as.character(dat4plot[,4]))
	
	p = ggplot(dat4plot, aes(x=cancer, y=value)) + geom_boxplot( aes(color = group)) + 
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjus=1), plot.title = element_text(hjust = 0.5) ) + ggtitle(drugs[kdrug]) + theme(legend.position="none") + 
		ylab("log(TMB)") + xlab("")
	print(p)	
}

p = ggplot(dat4plot, aes(x=cancer, y=value, fill=group)) + geom_boxplot( aes(color = group)) 
print(p)	
dev.off()
