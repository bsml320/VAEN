setwd("/path/to/GitHub/Figure/Figure5/GSE33072")

library("GEOquery")
library("limma")

gse <- getGEO("GSE33072", GSEMatrix = TRUE)

save(gse, file="GSE33072.RData")

exprs(gse[[1]]) -> datExpr0

########## section 1

library('preprocessCore')
normalize.quantiles(datExpr0) -> datExprLQ
dimnames(datExprLQ) = dimnames(datExpr0)

########## section 1 end


########## section 2

gpl = getGEO("GPL6244")
Table(gpl) -> anno

anno$gene_assignment -> ss
sapply(ss, function(u)trimws(strsplit(u, split="//")[[1]][2])) -> x

anno = anno[which(!is.na(x)), ]
anno$gene_assignment -> ss
sapply(ss, function(u)trimws(strsplit(u, split="//")[[1]][2])) -> x
names(x) = NULL

datExprLQ = datExprLQ[which(rownames(datExprLQ) %in% anno[,1]),]

match(rownames(datExprLQ), anno[,1]) -> ii
symbol = x[ii]

apply(datExprLQ, 2, function(u)tapply(u, symbol, median)) -> expr.mat

pData(gse[[1]]) -> pheno.anno

save(expr.mat, pheno.anno, file="GSE33072.RData")
