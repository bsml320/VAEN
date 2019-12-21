setwd("/work/Figures/Figure9/")

###################################################################################################################

drug.ccle = read.table(file="/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

cancer.drug.ccle = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	
	type.code = "01"
	if(cancer == "LAML"){ type.code = "03" }
	if(cancer == "SKCM"){ type.code = "06" }
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
	cancer.drug.ccle = rbind(cancer.drug.ccle, blca.ccle)
}
drug.ccle = cancer.drug.ccle


###################################################################################################################
stage.res = c()
cancer.type = unique(drug.ccle[,2])
for(cancer in cancer.type){
	clin.data = read.delim(paste("/work/data/TCGA/",cancer,"/",cancer,"_clinicalMatrix", sep=""),as.is=T)
	substr(clin.data[,1], 14, 15) -> type
	code = "01"
	if(cancer == "SKCM")code = "06"
	if(cancer == "LAML")code = "03"
	clin.data = clin.data[type==code,]
	
	stage.col.idx = 0
	if(is.element("pathologic_stage", colnames(clin.data))){
		stage.col.idx = which(colnames(clin.data) == "pathologic_stage")
	} else if(is.element("clinical_stage", colnames(clin.data))){
		stage.col.idx = which(colnames(clin.data) == "clinical_stage")
	} else {
		print(cancer)
		next
	}
	
	which(clin.data[, stage.col.idx] %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC", "Stage IA1", "Stage IA2", "Stage IB1", "Stage IB2") ) -> grp.I.ii
	which(clin.data[, stage.col.idx] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC", "Stage IIA1", "Stage IIA2") ) -> grp.II.ii
	which(clin.data[, stage.col.idx] %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ) -> grp.III.ii
	which(clin.data[, stage.col.idx] %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage X") ) -> grp.IV.ii

	which(drug.ccle[,1] %in% clin.data[grp.I.ii, 1]) -> ii.1
	which(drug.ccle[,1] %in% clin.data[grp.II.ii, 1]) -> ii.2
	which(drug.ccle[,1] %in% clin.data[grp.III.ii, 1]) -> ii.3
	which(drug.ccle[,1] %in% clin.data[grp.IV.ii, 1]) -> ii.4
	
	if(length(ii.4) < 5){
		print(cancer)
		next
	}
	
	pdf(paste(cancer, ".stage.pdf", sep=""))
	for(drug in colnames(drug.ccle)[c(-1,-2)]){
		wilcox.test(drug.ccle[c(ii.1, ii.2, ii.3), drug], drug.ccle[c(ii.4), drug])$p.value -> p
		boxplot(list(I = drug.ccle[ii.1, drug], II=drug.ccle[ii.2, drug], III=drug.ccle[ii.3, drug], IV=drug.ccle[ii.4, drug]), main=paste(drug, ", p = ", format(p, digits=3), sep=""))
		stage.res = rbind(stage.res, c(cancer, drug, length(ii.1), length(ii.2), length(ii.3), length(ii.4), p, ifelse(mean(drug.ccle[c(ii.1, ii.2, ii.3), drug]) > mean(drug.ccle[ii.4, drug]), 1, -1) ) )
	}
	dev.off()

}
write.table(stage.res, file="CCLE.stage.res.txt", row.names=F, col.names=F, quote=F, sep="\t")

library("reshape2")
dat = read.table("/work/Figures/Figure9/CCLE.stage.res.txt", as.is=T)
acast(dat, V1~V2, value.var="V7") -> test

dat = cbind(dat, log10(dat[,7]) * dat[,8])
colnames(dat)[9] = "V9"
acast(dat, V1~V2, value.var="V9") -> test


library("pheatmap")
test_label = test
temp = abs(test_label) > 1.3
test_label[temp] <- "*"
test_label[!temp] <- ""

pdf("9B.stage.pdf")
pheatmap(test, display_numbers = test_label,col = c(colorRampPalette(c("green", "white"))(floor(min(test)*(-10))), c(colorRampPalette(c("white", "red"))(floor(max(test)*(10))))), fontsize_number=20)
dev.off()



