#setwd("/path/to/VAEN/Figure/FigureS13")

################## prepare cancer mutations
load("../Figure7/MC3.RData")
cancer.types = names(cancer.mutation.list)

################################################################################################################

drug.ccle = read.table(file="../../result.EN/dr.GDSC/VAEN_GDSC.A.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

#########################################################################################################

cluster.drug.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	if(cancer == "BRCA" | cancer == "LGG" | cancer == "THCA")next
	mut.list = cancer.mutation.list[[ cancer ]]
	fixed.ss = mut.list[[ "ss"  ]]
	gene2WT.list = mut.list[[ "gene2WT.list" ]]
	new.mut.mat.3cols = read.table(paste("../Figure7/cluster/",cancer, ".dSNP.4cluster.MC3.txt", sep=""), as.is=T, header=T)
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	
	tapply(new.mut.mat.3cols$Tumor_Sample_Barcode3, new.mut.mat.3cols$Gene, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = new.mut.mat.3cols[which(new.mut.mat.3cols$Gene==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = gene2WT.list[[ symbol ]]
		as.numeric(names(table(gene.mutation[,4]))) -> mut.positions
		g = rep(1, length(mut.positions))
		if(length(mut.positions)>1){
			for(j in 2:length(mut.positions)){
				if(mut.positions[j] - mut.positions[j-1] <= 5){
					g[j] = g[j-1]
				} else {
					g[j] = g[j-1] + 1
				}
			}
		}
		names(g) = as.character(mut.positions)
		
		###################################################################################
		### X: mutation position
		cutoff = 10
		names( table(g) ) -> mut.clusters
		X.mut_pos.mat = matrix(NA, nrow=length(fixed.ss), ncol=length(mut.clusters))
		rownames(X.mut_pos.mat) = fixed.ss
		for(gg in 1:length(mut.clusters)){
			names(which(g==mut.clusters[gg])) -> actual.mut_pos
			which(as.character(gene.mutation$mut_pos) %in% actual.mut_pos) -> g.ss.idx
			X.mut_pos.mat[gene.mutation$Tumor_Sample_Barcode3[g.ss.idx], gg] = as.numeric(mut.clusters[gg])
			X.mut_pos.mat[WT.ss, gg] = 0
		}
		colnames(X.mut_pos.mat) = mut.clusters
		
		apply(X.mut_pos.mat, 2, function(u)sum(u!=0, na.rm=T)) -> count.each.cluster
		if(sum(count.each.cluster>=cutoff)!=0){
			for(kc in 1:ncol(X.mut_pos.mat)){
				u = X.mut_pos.mat[, kc]
				if(sum(u!=0, na.rm=T) < cutoff)next
				paste(names(g[which(g==colnames(X.mut_pos.mat)[kc])]), collapse=":") -> clust.name
				
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					wilcox.test(v[  which(u==colnames(X.mut_pos.mat)[kc])  ], v[which(u==0)])$p.value 
				})-> ps.twosided
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					mean(v[u!=0], na.rm=T)-mean(v[which(u==0)])
				})-> fcs
				cluster.drug.mat = rbind(cluster.drug.mat, cbind(cancer, symbol, gene, colnames(blca.ccle)[3:ncol(blca.ccle)], names(ps.twosided), ps.twosided, fcs, nMut = sum(u!=0, na.rm=T), clust.name, nMut.prop = sum(u!=0, na.rm=T)/nrow(blca.ccle) ))
			}
		}
	}
	cat(cancer, ".", sep="")
}

####################################################################
cancer = "LGG"
mut.list = cancer.mutation.list[[ cancer ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]
new.mut.mat.3cols = read.table(paste("../Figure7/cluster/", cancer, ".dSNP.4cluster.MC3.txt", sep=""), as.is=T, header=T)

subtype = read.delim("../../DATA/LGG.subtype.TCGA.txt", as.is=T, header=T)
ss = substr(fixed.ss, 1, 12)
subtype[which(subtype[,9]=="coc1"),1] -> g1
subtype[which(subtype[,9]=="coc2"),1] -> g2
subtype[which(subtype[,9]=="coc3"),1] -> g3
subtype.list = list()
subtype.list[["LGG-coc1"]] = g1
subtype.list[["LGG-coc2"]] = g2
subtype.list[["LGG-coc3"]] = g3
	
for(k in 1:length(subtype.list)){
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	substr(new.mut.mat.3cols[,2], 1, 12) -> new.mut.mat.3cols.ss
	subtype.mut.mat = new.mut.mat.3cols[new.mut.mat.3cols.ss %in% subtype.ss, ]
	fixed.ss = unique(subtype.mut.mat[,2])
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	
	tapply(subtype.mut.mat$Tumor_Sample_Barcode3, subtype.mut.mat$Gene, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = subtype.mut.mat[which(subtype.mut.mat$Gene==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = intersect(gene2WT.list[[ symbol ]], fixed.ss)
		as.numeric(names(table(gene.mutation[,4]))) -> mut.positions
		g = rep(1, length(mut.positions))
		if(length(mut.positions)>1){
			for(j in 2:length(mut.positions)){
				if(mut.positions[j] - mut.positions[j-1] <= 5){
					g[j] = g[j-1]
				} else {
					g[j] = g[j-1] + 1
				}
			}
		}
		names(g) = as.character(mut.positions)
		
		###################################################################################
		### X: mutation position
		cutoff = 10
		names( table(g) ) -> mut.clusters
		X.mut_pos.mat = matrix(NA, nrow=length(fixed.ss), ncol=length(mut.clusters))
		rownames(X.mut_pos.mat) = fixed.ss
		for(gg in 1:length(mut.clusters)){
			names(which(g==mut.clusters[gg])) -> actual.mut_pos
			which(as.character(gene.mutation$mut_pos) %in% actual.mut_pos) -> g.ss.idx
			X.mut_pos.mat[gene.mutation$Tumor_Sample_Barcode3[g.ss.idx], gg] = as.numeric(mut.clusters[gg])
			X.mut_pos.mat[WT.ss, gg] = 0
		}
		colnames(X.mut_pos.mat) = mut.clusters
		
		apply(X.mut_pos.mat, 2, function(u)sum(u!=0, na.rm=T)) -> count.each.cluster
		if(sum(count.each.cluster>=cutoff)!=0){
			for(kc in 1:ncol(X.mut_pos.mat)){
				u = X.mut_pos.mat[, kc]
				if(sum(u!=0, na.rm=T) < cutoff)next
				paste(names(g[which(g==colnames(X.mut_pos.mat)[kc])]), collapse=":") -> clust.name
				
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					wilcox.test(v[  which(u==colnames(X.mut_pos.mat)[kc])  ], v[which(u==0)])$p.value 
				})-> ps.twosided
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					mean(v[u!=0], na.rm=T)-mean(v[which(u==0)])
				})-> fcs
				cluster.drug.mat = rbind(cluster.drug.mat, cbind(cancer, symbol, gene, colnames(blca.ccle)[3:ncol(blca.ccle)], names(ps.twosided), ps.twosided, fcs, nMut = sum(u!=0, na.rm=T), clust.name, nMut.prop = sum(u!=0, na.rm=T)/nrow(blca.ccle) ))
			}
		}
	}
	cat(cancer, ".", sep="")
}

####################################################################
cancer = "BRCA"
mut.list = cancer.mutation.list[[ cancer ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]
new.mut.mat.3cols = read.table(paste("../Figure7/cluster/", cancer, ".dSNP.4cluster.MC3.txt", sep=""), as.is=T, header=T)

subtype = read.delim("../../DATA/BRCA.subtype.TCGA.txt", as.is=T, header=T, sep="\t", skip=1)
ss = substr(fixed.ss, 1, 12)
subtype[which(subtype$PAM50.mRNA=="Basal-like"),1] -> g1
subtype[which(subtype$PAM50.mRNA=="HER2-enriched"),1] -> g2
subtype[which(subtype$PAM50.mRNA=="Luminal A"),1] -> g3
subtype[which(subtype$PAM50.mRNA=="Luminal B"),1] -> g4
subtype.list = list()
subtype.list[["BRCA-Basal"]] = g1
subtype.list[["BRCA-HER2"]] = g2
subtype.list[["BRCA-LumA"]] = g3
subtype.list[["BRCA-LumB"]] = g4
	
for(k in 1:length(subtype.list)){
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	substr(new.mut.mat.3cols[,2], 1, 12) -> new.mut.mat.3cols.ss
	subtype.mut.mat = new.mut.mat.3cols[new.mut.mat.3cols.ss %in% subtype.ss, ]
	fixed.ss = unique(subtype.mut.mat[,2])
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	
	tapply(subtype.mut.mat$Tumor_Sample_Barcode3, subtype.mut.mat$Gene, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = subtype.mut.mat[which(subtype.mut.mat$Gene==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = intersect(gene2WT.list[[ symbol ]], fixed.ss)
		as.numeric(names(table(gene.mutation[,4]))) -> mut.positions
		g = rep(1, length(mut.positions))
		if(length(mut.positions)>1){
			for(j in 2:length(mut.positions)){
				if(mut.positions[j] - mut.positions[j-1] <= 5){
					g[j] = g[j-1]
				} else {
					g[j] = g[j-1] + 1
				}
			}
		}
		names(g) = as.character(mut.positions)
		
		###################################################################################
		### X: mutation position
		cutoff = 10
		names( table(g) ) -> mut.clusters
		X.mut_pos.mat = matrix(NA, nrow=length(fixed.ss), ncol=length(mut.clusters))
		rownames(X.mut_pos.mat) = fixed.ss
		for(gg in 1:length(mut.clusters)){
			names(which(g==mut.clusters[gg])) -> actual.mut_pos
			which(as.character(gene.mutation$mut_pos) %in% actual.mut_pos) -> g.ss.idx
			X.mut_pos.mat[gene.mutation$Tumor_Sample_Barcode3[g.ss.idx], gg] = as.numeric(mut.clusters[gg])
			X.mut_pos.mat[WT.ss, gg] = 0
		}
		colnames(X.mut_pos.mat) = mut.clusters
		
		apply(X.mut_pos.mat, 2, function(u)sum(u!=0, na.rm=T)) -> count.each.cluster
		if(sum(count.each.cluster>=cutoff)!=0){
			for(kc in 1:ncol(X.mut_pos.mat)){
				u = X.mut_pos.mat[, kc]
				if(sum(u!=0, na.rm=T) < cutoff)next
				paste(names(g[which(g==colnames(X.mut_pos.mat)[kc])]), collapse=":") -> clust.name
				
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					wilcox.test(v[  which(u==colnames(X.mut_pos.mat)[kc])  ], v[which(u==0)])$p.value 
				})-> ps.twosided
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					mean(v[u!=0], na.rm=T)-mean(v[which(u==0)])
				})-> fcs
				cluster.drug.mat = rbind(cluster.drug.mat, cbind(cancer, symbol, gene, colnames(blca.ccle)[3:ncol(blca.ccle)], names(ps.twosided), ps.twosided, fcs, nMut = sum(u!=0, na.rm=T), clust.name, nMut.prop = sum(u!=0, na.rm=T)/nrow(blca.ccle) ))
			}
		}
	}
	cat(cancer, ".", sep="")
}

####################################################################
### THCA
cancer = "THCA"
mut.list = cancer.mutation.list[[ cancer ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]
new.mut.mat.3cols = read.table(paste("../Figure7/cluster/", cancer, ".dSNP.4cluster.MC3.txt", sep=""), as.is=T, header=T)

subtype = read.delim("../../DATA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
ss = substr(subtype[,1], 1, 15)
ss[which(subtype[,2]=="Braf-like")] -> g1
ss[which(subtype[,2]=="Ras-like")] -> g2
subtype.list = list()
subtype.list[["THCA-Braf-like"]] = g1
subtype.list[["THCA-Ras-like"]] = g2
	
for(k in 1:length(subtype.list)){
	cancer = names(subtype.list)[k]
	subtype.ss = subtype.list[[k]]
	substr(new.mut.mat.3cols[,2], 1, 12) -> new.mut.mat.3cols.ss
	subtype.mut.mat = new.mut.mat.3cols[new.mut.mat.3cols.ss %in% subtype.ss, ]
	fixed.ss = unique(subtype.mut.mat[,2])
	
	blca.ccle = drug.ccle[match(fixed.ss, drug.ccle[,1]),]
	
	tapply(subtype.mut.mat$Tumor_Sample_Barcode3, subtype.mut.mat$Gene, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 10)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = subtype.mut.mat[which(subtype.mut.mat$Gene==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = intersect(gene2WT.list[[ symbol ]], fixed.ss)
		as.numeric(names(table(gene.mutation[,4]))) -> mut.positions
		g = rep(1, length(mut.positions))
		if(length(mut.positions)>1){
			for(j in 2:length(mut.positions)){
				if(mut.positions[j] - mut.positions[j-1] <= 5){
					g[j] = g[j-1]
				} else {
					g[j] = g[j-1] + 1
				}
			}
		}
		names(g) = as.character(mut.positions)
		
		###################################################################################
		### X: mutation position
		cutoff = 10
		names( table(g) ) -> mut.clusters
		X.mut_pos.mat = matrix(NA, nrow=length(fixed.ss), ncol=length(mut.clusters))
		rownames(X.mut_pos.mat) = fixed.ss
		for(gg in 1:length(mut.clusters)){
			names(which(g==mut.clusters[gg])) -> actual.mut_pos
			which(as.character(gene.mutation$mut_pos) %in% actual.mut_pos) -> g.ss.idx
			X.mut_pos.mat[gene.mutation$Tumor_Sample_Barcode3[g.ss.idx], gg] = as.numeric(mut.clusters[gg])
			X.mut_pos.mat[WT.ss, gg] = 0
		}
		colnames(X.mut_pos.mat) = mut.clusters
		
		apply(X.mut_pos.mat, 2, function(u)sum(u!=0, na.rm=T)) -> count.each.cluster
		if(sum(count.each.cluster>=cutoff)!=0){
			for(kc in 1:ncol(X.mut_pos.mat)){
				u = X.mut_pos.mat[, kc]
				if(sum(u!=0, na.rm=T) < cutoff)next
				paste(names(g[which(g==colnames(X.mut_pos.mat)[kc])]), collapse=":") -> clust.name
				
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					wilcox.test(v[  which(u==colnames(X.mut_pos.mat)[kc])  ], v[which(u==0)])$p.value 
				})-> ps.twosided
				apply(blca.ccle[, 3:ncol(blca.ccle)], 2, function(v){
					mean(v[u!=0], na.rm=T)-mean(v[which(u==0)])
				})-> fcs
				cluster.drug.mat = rbind(cluster.drug.mat, cbind(cancer, symbol, gene, colnames(blca.ccle)[3:ncol(blca.ccle)], names(ps.twosided), ps.twosided, fcs, nMut = sum(u!=0, na.rm=T), clust.name, nMut.prop = sum(u!=0, na.rm=T)/nrow(blca.ccle) ))
			}
		}
		###################################################################################
	}
	cat(cancer, ".", sep="")
}

write.table(cluster.drug.mat, file=paste("tmp", sep=""), quote=F, sep="\t", row.names=F, col.names=F )
dSNP.cluster.mat = read.table("tmp", header=F, as.is=T, sep="\t")
file.remove("tmp")
dSNP.cluster.mat = dSNP.cluster.mat[, -5]

###### *****
drugs = unique(dSNP.cluster.mat[,4])
dSNP.cluster.mat$adjp = 1
for(k in 1:length(drugs)){
	which(dSNP.cluster.mat[,4] == drugs[k]) -> ii
	p.adjust(dSNP.cluster.mat[ii,5], method="BH") -> adjp
	dSNP.cluster.mat[ii, "adjp"] = adjp
}
###### *****

which(dSNP.cluster.mat$adjp < 0.05 & dSNP.cluster.mat[,6] > 0) -> ii.1
which(dSNP.cluster.mat$adjp < 0.05 & dSNP.cluster.mat[,6] < 0) -> ii.2
cc = rep(0, nrow(dSNP.cluster.mat)); cc[ii.1] = 2; cc[ii.2] = 1
dSNP.cluster.mat$color = as.factor(cc)

write.table(dSNP.cluster.mat, file="FigureS13.B.txt", row.names=F, quote=F, sep="\t")
