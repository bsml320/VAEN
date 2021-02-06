##############################################################################################################
### For Figure 3A
setwd("/path/to/VAEN/Figure/Figure3")

### load dr.ccle.models
load("../../result.EN/dr.CCLE/dr.CCLE.A.models.RData")
drugs = names(dr.ccle.models)

TCGA.pred = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")

#####################
dr.ccle.mat = c()
for(kdrug in 1:length(drugs)){
	drug = drugs[kdrug]
	
	dr.ccle.models[[ drug ]] -> res.list
	Ys = res.list$Ys
	
	which(Ys[,1]!=-9) -> ii
	Ys = Ys[ii, ]
	dr.ccle.mat = rbind(dr.ccle.mat, cbind(drug = drug, sample = rownames(Ys), ActArea = Ys[,1], grp="CCLE Observed DR") )
	dr.ccle.mat = rbind(dr.ccle.mat, cbind(drug = drug, sample = rownames(Ys), ActArea = Ys[,2], grp="CCLE Predicted DR") )
	
	TCGA.drug.name = gsub("-", ".", drug)
	if(drug == "17-AAG")TCGA.drug.name = "X17.AAG"
	cur.drug.TCGA.pred = TCGA.pred[, TCGA.drug.name]
	dr.ccle.mat = rbind(dr.ccle.mat, cbind(drug = drug, sample = gsub("-",".",TCGA.pred[,1]), ActArea = cur.drug.TCGA.pred, grp="TCGA Predicted DR") )
}

#########################

new.mat = as.data.frame(dr.ccle.mat)
rownames(new.mat) = NULL
new.mat[,3] = as.numeric(as.character(new.mat[,3]))
write.table(new.mat, file="Figure3A.txt", row.names=F, quote=F, sep="\t")

pdf("3A.pdf", width=10, height=5)
g = ggplot(aes(y = ActArea, x = drug, fill = grp), data = new.mat) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(0.9, 0.9)) + 
  xlab("") + ylab("Drug response (ActArea)") +
  guides(fill=guide_legend(title=""))
print(g)
dev.off()


