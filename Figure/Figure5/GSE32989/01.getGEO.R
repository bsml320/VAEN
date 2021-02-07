setwd("/path/to/GitHub/result.EN/dr.CCLE/FIGURE/FIGURE5/GSE32989")
library("GEOquery")
library("limma")

gse <- getGEO("GSE32989", GSEMatrix = TRUE)

save(gse, file="GSE32989.RData")

exprs(gse[[1]]) -> datExpr0

########## section 1

library('preprocessCore')
normalize.quantiles(datExpr0) -> datExprLQ
dimnames(datExprLQ) = dimnames(datExpr0)

########## section 1 end


########## section 2

gpl = getGEO("GPL13376")
Table(gpl) -> anno

match(rownames(datExprLQ), anno[,1]) -> ii
symbol = anno$ILMN_Gene[ii]

apply(datExprLQ, 2, function(u)tapply(u, symbol, median)) -> expr.mat

pData(gse[[1]]) -> pheno.anno

save(expr.mat, pheno.anno, file="GSE32989.RData")
