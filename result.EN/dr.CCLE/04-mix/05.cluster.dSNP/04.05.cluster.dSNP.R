setwd("/work/result.EN/dr.CCLE/04-mix/05.cluster.dSNP/")

################## prepare cancer mutations
load("/work/result.EN/dr.CCLE/04-mix/01.dSNP/cancer.mutation.list.dSNP.WT.RData")
cancer.types = names(cancer.mutation.list)

for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	mut.list = cancer.mutation.list[[ cancer ]]
	new.mut.mat.3cols = mut.list[[ "mut" ]]
	
	apply(new.mut.mat.3cols, 1, function(u){
		strsplit(u[3], split=",")[[1]] -> v
		new.record = c()
		for(k in 1:length(v)){
			strsplit(v[k], split=":")[[1]] -> v1
			if(length(v1) == 5){
				new.record = rbind(new.record, c(u[1], u[2], unlist(v1)) )
			} else {
				print(u[1])
			}
		}
		new.record
	}) -> mut.split.list
	do.call(rbind, mut.split.list) -> x
	
	as.numeric(unlist( lapply(x[,7], function(u){ 
		if(grepl("_", u)){
			substr(u, 3, regexec("_", u)[[1]][1]-1  ) -> u1
			regexec("[0-9]+", u1) -> m
			regmatches(u1,m) -> a
		} else {
			regexec("[0-9]+", u) -> m
			regmatches(u,m) -> a
		} } ) )) -> mut_pos
	new.mut.mat.8cols = cbind(x, mut_pos)
	colnames(new.mut.mat.8cols) = c("Gene", "Tumor_Sample_Barcode3", "Gene2", "mRNA_ID", "Exon", "nucleotide", "AAchange", "mut_pos")
	write.table(new.mut.mat.8cols, file=paste(cancer, ".dSNP.4cluster.txt", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
	cat(cancer,".",sep="")
}

################################################################################################################

drug.ccle = read.table(file="/work/result.EN/dr.CCLE/01/MIX-F1-W5-PCC/MIX-F1-W5-PCC.avgtop10.pred_TCGA.txt", header=T, as.is=T, sep="\t")
colnames(drug.ccle)[3:ncol(drug.ccle)] -> drugs
cancer.types = unique(drug.ccle[,2])
sample.type = substr(drug.ccle[,1], 14, 15)
ss = gsub("\\.", "-", drug.ccle[,1])
drug.ccle[,1] = ss

###################################

cluster.drug.mat = c()
for(ct in 1:length(cancer.types)){
	cancer = cancer.types[ct]
	mut.list = cancer.mutation.list[[ cancer ]]
	fixed.ss = mut.list[[ "ss"  ]]
	gene2WT.list = mut.list[[ "gene2WT.list" ]]
	new.mut.mat.3cols = read.table(paste(cancer, ".dSNP.4cluster.txt", sep=""), as.is=T, header=T)
	
	type.code = "01"
	if(cancer == "LAML"){type.code = "03"}
	if(cancer == "SKCM"){type.code = "06"}
	
	blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type==type.code), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	
	tapply(new.mut.mat.3cols$Tumor_Sample_Barcode3, new.mut.mat.3cols$mRNA_ID, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = new.mut.mat.3cols[which(new.mut.mat.3cols$mRNA_ID==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = gene2WT.list[[ symbol ]]
		as.numeric(names(table(gene.mutation[,8]))) -> mut.positions
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

####################################################################
cancer = "LGG"
mut.list = cancer.mutation.list[[ cancer ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]
new.mut.mat.3cols = read.table(paste(cancer, ".dSNP.4cluster.txt", sep=""), as.is=T, header=T)

type.code = "01"
blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle

subtype = read.delim("/work/data/TCGA/LGG.subtype.TCGA.txt", as.is=T, header=T)
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
	
	type.code = "01"
	blca.ccle = drug.ccle[which(drug.ccle[,2] == "LGG" & sample.type==type.code), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	
	tapply(subtype.mut.mat$Tumor_Sample_Barcode3, subtype.mut.mat$mRNA_ID, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = subtype.mut.mat[which(subtype.mut.mat$mRNA_ID==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = intersect(gene2WT.list[[ symbol ]], fixed.ss)
		as.numeric(names(table(gene.mutation[,8]))) -> mut.positions
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

####################################################################
cancer = "BRCA"
mut.list = cancer.mutation.list[[ cancer ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]
new.mut.mat.3cols = read.table(paste(cancer, ".dSNP.4cluster.txt", sep=""), as.is=T, header=T)

type.code = "01"
blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle

subtype = read.delim("/work/data/TCGA/BRCA.subtype.TCGA.txt", as.is=T, header=T, sep="\t", skip=1)
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
	
	type.code = "01"
	blca.ccle = drug.ccle[which(drug.ccle[,2] == "BRCA" & sample.type==type.code), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	
	tapply(subtype.mut.mat$Tumor_Sample_Barcode3, subtype.mut.mat$mRNA_ID, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = subtype.mut.mat[which(subtype.mut.mat$mRNA_ID==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = intersect(gene2WT.list[[ symbol ]], fixed.ss)
		as.numeric(names(table(gene.mutation[,8]))) -> mut.positions
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

####################################################################
### THCA
cancer = "THCA"
mut.list = cancer.mutation.list[[ cancer ]]
fixed.ss = mut.list[[ "ss"  ]]
gene2WT.list = mut.list[[ "gene2WT.list" ]]
new.mut.mat.3cols = read.table(paste(cancer, ".dSNP.4cluster.txt", sep=""), as.is=T, header=T)

type.code = "01"
blca.ccle = drug.ccle[which(drug.ccle[,2] == cancer & sample.type == type.code), ]
blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
apply(blca.ccle[, 3:ncol(blca.ccle)], 2, as.numeric) -> new.ccle

subtype = read.delim("/work/data/TCGA/THCA.BRS.txt", as.is=T, header=F, sep="\t")
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
	
	type.code = "01"
	blca.ccle = drug.ccle[which(drug.ccle[,2] == "BRCA" & sample.type==type.code), ]
	blca.ccle = blca.ccle[match(fixed.ss, blca.ccle[,1]),]
	
	tapply(subtype.mut.mat$Tumor_Sample_Barcode3, subtype.mut.mat$mRNA_ID, function(u)length(unique(u))) -> gene2ss
	candidate.genes = sort(names(gene2ss[which(gene2ss >= 5)]))
	if(length(candidate.genes) < 1)next
	
	for(kgene in 1:length(candidate.genes)){
		gene = candidate.genes[kgene]
		gene.mutation = subtype.mut.mat[which(subtype.mut.mat$mRNA_ID==gene), ]
		symbol = gene.mutation[1,1]
		WT.ss = intersect(gene2WT.list[[ symbol ]], fixed.ss)
		as.numeric(names(table(gene.mutation[,8]))) -> mut.positions
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

write.table(cluster.drug.mat, file=paste("04.05.cluster.dSNP.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=F )

