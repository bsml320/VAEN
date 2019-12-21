setwd("/work/Figures/Figure9")

###################################################################################################################

drug.ccle = read.table(file="/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
drugs = colnames(drug.ccle)[c(-1, -2)]
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)

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

dat4boxplot.dat = c()
cancer.type = unique(drug.ccle[,2])
drug = "Irinotecan"
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
	
	dat4boxplot.dat = rbind(dat4boxplot.dat, cbind(cancer=cancer, drug=drug, group="Stage IV", pred.dr = drug.ccle[ii.4, drug]), cbind(cancer=cancer, drug=drug, group="Stage I, II, III", pred.dr = drug.ccle[c(ii.1, ii.2, ii.3), drug])  )
}

#
#[1] "GBM"
#[1] "LAML"
#[1] "LGG"
#[1] "PAAD"
#[1] "PCPG"
#[1] "PRAD"
#[1] "SARC"
#[1] "TGCT"
#[1] "THYM"
#[1] "UVM"
#


dat4boxplot.dat = as.data.frame(dat4boxplot.dat)
dat4boxplot.dat[,4] = as.numeric(as.character(dat4boxplot.dat[,4]))
dat4boxplot.dat[,1] <- factor(dat4boxplot.dat[,1], levels = c("ACC", "KICH", "KIRC", "KIRP", "CESC", "UCEC", "UCS", "BLCA", "BRCA", "COAD", "CHOL", "DLBC", "ESCA", "HNSC", "LIHC", "LUAD", "LUSC", "MESO", "OV", "READ", "SKCM", "STAD", "THCA"))


p1 = ggplot(dat4boxplot.dat,aes(x=cancer ,y=pred.dr,color=group,fill=group))  + 
	  geom_point(position=position_jitterdodge(dodge.width=1), size = 1, stroke = 0, shape = 16) + 
	  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9), lwd=0.4, fatten=0.8) + 
	  xlab("") + ylab("Imputed Drug Response") + ggtitle("Irinotecan") + 
	  theme(plot.title = element_text(hjust = 0.5, size=8), legend.position = c(0.9, 0.9), axis.text.x = element_text(angle = 30, hjust = 1, size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6) )
      
pdf("9A.Irinotecan.pdf", width=5, height=2.2)
print(p1)	  
dev.off()

############################################################################################

dat4boxplot.dat = c()
cancer.type = unique(drug.ccle[,2])
drug = "Topotecan"
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
	
	dat4boxplot.dat = rbind(dat4boxplot.dat, cbind(cancer=cancer, drug=drug, group="Stage IV", pred.dr = drug.ccle[ii.4, drug]), cbind(cancer=cancer, drug=drug, group="Stage I, II, III", pred.dr = drug.ccle[c(ii.1, ii.2, ii.3), drug])  )
}


dat4boxplot.dat = as.data.frame(dat4boxplot.dat)
dat4boxplot.dat[,4] = as.numeric(as.character(dat4boxplot.dat[,4]))
dat4boxplot.dat[,1] <- factor(dat4boxplot.dat[,1], levels = c("ACC", "KICH", "KIRC", "KIRP", "CESC", "UCEC", "UCS", "BLCA", "BRCA", "COAD", "CHOL", "DLBC", "ESCA", "HNSC", "LIHC", "LUAD", "LUSC", "MESO", "OV", "READ", "SKCM", "STAD", "THCA"))

p1 = ggplot(dat4boxplot.dat,aes(x=cancer ,y=pred.dr,color=group,fill=group))  + 
	  geom_point(position=position_jitterdodge(dodge.width=1), size = 1, stroke = 0, shape = 16) + 
	  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9), lwd=0.4, fatten=0.8) + 
	  xlab("") + ylab("Imputed Drug Response") + ggtitle("Topotecan") + 
	  theme(plot.title = element_text(hjust = 0.5, size=8), legend.position = c(0.9, 0.9), axis.text.x = element_text(angle = 30, hjust = 1, size=6), axis.text.y = element_text(size=6), axis.title.y=element_text(size=6), legend.text=element_text(size=6), legend.title=element_text(size=6) )
      
pdf("9A.Topotecan.pdf", width=5, height=2.2)
print(p1)	  
dev.off()

#################################################################################################


pmat = c()
cancer.type = unique(drug.ccle[,2])
for(cancer in cancer.type){
	clin.data = read.delim(paste("/data/mshao/TCGA/",cancer,"/",cancer,"_clinicalMatrix", sep=""),as.is=T)
	substr(clin.data[,1], 14, 15) -> type
	code = "01"
	if(cancer == "SKCM")code = "06"
	if(cancer == "LAML")code = "03"
	clin.data = clin.data[type==code,]
	clin.data[, 1] = substr(clin.data[,1], 1, 12)
	
	stage.col.idx = 0
	if(is.element("pathologic_stage", colnames(clin.data))){
		stage.col.idx = which(colnames(clin.data) == "pathologic_stage")
	} else if(is.element("clinical_stage", colnames(clin.data))){
		stage.col.idx = which(colnames(clin.data) == "clinical_stage")
	} else {
		print(cancer)
		next
	}
	
	fn = paste("/data1_2/jiap/projects/18-CCLE-VAE/new/V9/ssGSEA/",cancer,".ssgsea.result.RData", sep="")
	if(!file.exists(fn))next
	load(fn)
	ssgsea.type = sapply(colnames(es.mat), function(u)substr(u, 14, 15)) -> type
	es.mat = es.mat[, which(ssgsea.type == type.code)]
	ssgsea.ss = gsub("\\.","-",substr(colnames(es.mat), 1, 12))
	colnames(es.mat) = ssgsea.ss
	
	
	which(clin.data[, stage.col.idx] %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC", "Stage IA1", "Stage IA2", "Stage IB1", "Stage IB2") ) -> grp.I.ii
	which(clin.data[, stage.col.idx] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC", "Stage IIA1", "Stage IIA2") ) -> grp.II.ii
	which(clin.data[, stage.col.idx] %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ) -> grp.III.ii
	which(clin.data[, stage.col.idx] %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage X") ) -> grp.IV.ii

	which( colnames(es.mat) %in% clin.data[grp.I.ii, 1]) -> ii.1
	which( colnames(es.mat) %in% clin.data[grp.II.ii, 1]) -> ii.2
	which( colnames(es.mat) %in% clin.data[grp.III.ii, 1]) -> ii.3
	which( colnames(es.mat) %in% clin.data[grp.IV.ii, 1]) -> ii.4
	
	if(length(ii.4) < 5){
		print(cancer)
		next
	}
	
	allp = apply(es.mat, 1, function(u)t.test(u[ii.4], u[c(ii.1, ii.2, ii.3)], alternative="greater")$p.value)
	pmat = rbind(pmat, c(cancer, allp))
}
cbind(pmat[,1], apply(pmat, 1, function(u)min(as.numeric(u[-1]))), apply(pmat, 1, function(u)colnames(pmat)[which.min(as.numeric(u[-1]))+1]))
 
res.list = list()
for(k in 1:nrow(pmat)){
	u = as.numeric(pmat[k, -1])
	names(u) = colnames(pmat)[-1]
	u = sort(u)
	res.list[[pmat[k, 1]]] = u[1:10]
}

##############################

pdf("9C.barplot.Irinotecan.pdf")
par(mfrow=c(4,1))
data.list = list()
for(cancer in cancer.type){
	clin.data = read.delim(paste("/data/mshao/TCGA/",cancer,"/",cancer,"_clinicalMatrix", sep=""),as.is=T)
	substr(clin.data[,1], 14, 15) -> type
	code = "01"
	if(cancer == "SKCM")code = "06"
	if(cancer == "LAML")code = "03"
	clin.data = clin.data[type==code,]
	clin.data[, 1] = substr(clin.data[,1], 1, 12)
	
	stage.col.idx = 0
	if(is.element("pathologic_stage", colnames(clin.data))){
		stage.col.idx = which(colnames(clin.data) == "pathologic_stage")
	} else if(is.element("clinical_stage", colnames(clin.data))){
		stage.col.idx = which(colnames(clin.data) == "clinical_stage")
	} else {
		print(cancer)
		next
	}
	
	fn = paste("/work//data/ssGSEA/",cancer,".ssgsea.result.RData", sep="")
	if(!file.exists(fn))next
	load(fn)
	ssgsea.type = sapply(colnames(es.mat), function(u)substr(u, 14, 15)) -> type
	es.mat = es.mat[, which(ssgsea.type == type.code)]
	ssgsea.ss = gsub("\\.","-",substr(colnames(es.mat), 1, 12))
	colnames(es.mat) = ssgsea.ss
	
	
	which(clin.data[, stage.col.idx] %in% c("Stage I", "Stage IA", "Stage IB", "Stage IC", "Stage IA1", "Stage IA2", "Stage IB1", "Stage IB2") ) -> grp.I.ii
	which(clin.data[, stage.col.idx] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC", "Stage IIA1", "Stage IIA2") ) -> grp.II.ii
	which(clin.data[, stage.col.idx] %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ) -> grp.III.ii
	which(clin.data[, stage.col.idx] %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage X") ) -> grp.IV.ii

	which( colnames(es.mat) %in% clin.data[grp.I.ii, 1]) -> ii.1
	which( colnames(es.mat) %in% clin.data[grp.II.ii, 1]) -> ii.2
	which( colnames(es.mat) %in% clin.data[grp.III.ii, 1]) -> ii.3
	which( colnames(es.mat) %in% clin.data[grp.IV.ii, 1]) -> ii.4
	
	if(length(ii.4) < 5){
		print(cancer)
		next
	}
	pathway = es.mat["REACTOME_G0_AND_EARLY_G1", ]
	drug.ccle[match(names(pathway), substr(drug.ccle[,1], 1, 12)), "Irinotecan"] -> drugr
	names(drugr) = names(pathway)
	drugr = sort(drugr)
	pathway = pathway[names(drugr)]
	
	st = rep(0, length(pathway))
	st[which(names(drugr) %in% clin.data[grp.I.ii, 1]) ] = 1
	st[which(names(drugr) %in% clin.data[grp.II.ii, 1]) ] = 2
	st[which(names(drugr) %in% clin.data[grp.III.ii, 1]) ] = 3
	st[which(names(drugr) %in% clin.data[grp.IV.ii, 1]) ] = 4
	
	dat.mat = cbind(pathway, drugr, st)
	dat.mat = dat.mat[order(st, drugr),]
	pathway = dat.mat[,1]; drugr = dat.mat[,2]; st = dat.mat[,3]
	
	data.list[[cancer]] = dat.mat
	
	plot(0, ylim = c(0,0.7), xlim=c(-1, length(pathway)), col="white", main=cancer, xlab="", ylab="")
	cc = colorRampPalette(c("darkgreen", "white", "red"))(length(drugr))
	for(k in 1:length(drugr)){
		rect(k-0.3, 0.3, k+0.3, 0.5, col=cc[rank(drugr)[k] ], border=cc[rank(drugr)[k] ])
	}
	
	cc = colorRampPalette(c("white", "blue"))(length(pathway))
	for(k in 1:length(pathway)){
		rect(k-0.3, 0, k+0.3, 0.2, col=cc[rank(pathway)[k] ], border=cc[rank(pathway)[k] ])
	}
	
	drugr2 = c()
	for(k in 1:4){
		which(st == k) -> ii
		mean(drugr[ii]) -> a
		drugr2 = c(drugr2, a)
	}
	drugr2 = c(drugr2, drugr)
	cc2 = colorRampPalette(c("darkgreen", "white", "red"))(length(drugr2))
	
	cc = c("darkgreen", "lightgreen", "pink", "red")
	for(k in 1:4){
		which(st == k) -> ii
		rect(min(ii)-0.3, 0.52, max(ii)+0.3, 0.62, col=cc2[rank(drugr2)[k]], border="black")
		text(median(ii), 0.65, paste("Stage ", k, sep="") )
	}
	which(st == 0) -> ii
	if(length(ii) > 0)rect(min(ii)-0.3, 0.52, max(ii)+0.3, 0.62, col="grey", border="black")
	
	text(0, 0.1, "Pathway\nActivity")
	text(0, 0.4, "Irinotecan")
	
}
dev.off()

