setwd("/work/result.EN/dr.CCLE/04-mix/09.Expr/DAVID")
library("GOstats")
library("RDAVIDWebService")
library("org.Hs.eg.db")
library("clusterProfiler")


outputfolder <- ("/work/result.EN/dr.CCLE/04-mix/09.Expr/DAVID/DAVID_enrichment_plot_GOBP/")

#----------------------------------------------------------------------------DAVID 6.8 enrichment

files = dir("W5_Q95/")
files = files[grep("\\.txt", files)]

for (kf in 1:length(files)) {
	data = read.delim(paste("W5_Q95/",files[kf],sep=""), head = F, as.is=T)
	output = sub("\\.txt", "", files[kf]  )
	
	shared_EI = bitr(data$V1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

	david<-DAVIDWebService(email="1430803743@qq.com", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

	uplist=addList(david, shared_EI[,2], idType="ENTREZ_GENE_ID", listName="shared_EI", listType="Gene")

	categary=c("GOTERM_BP_DIRECT")
	setAnnotationCategories(david, categary)
	FuncAnnotChart <- getFunctionalAnnotationChart(david)
	
	write.table(FuncAnnotChart[,-6], file = paste(outputfolder, output,".txt",sep=""), row.names = F, quote=FALSE, sep="\t")
}
