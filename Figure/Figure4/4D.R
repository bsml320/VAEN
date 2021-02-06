setwd("/path/to/VAEN/Figure/Figure4")

give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

###################################################################################################
drug.ccle = read.table("../../result.EN/dr.CCLE/VAEN_CCLE.A.pred_TCGA.txt", header=T, as.is=T)
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

###################################################################################################

response = read.delim("../../DATA/response/drug_response.txt", as.is=T)

#############################
response = response[which(response$drug.name == "Paclitaxel"), ]

match(response[,2], substr(drug.ccle[,1], 1, 12)) -> ii
cbind(response, drug.ccle[ii, ] ) -> new2
new2 = new2[!is.na(ii), ]
dim(new2)

tapply(new2$Paclitaxel, new2$response, mean)
grep("Disease", new2$response) -> ii
drug = "Paclitaxel"
p = t.test(new2[ii, drug], new2[-ii, drug])$p.value


dat4plot = data.frame(new2[, c("Paclitaxel", "response")])
dat4plot[,1] = as.numeric(as.character(dat4plot[,1]))
dat4plot[,2] <- factor(dat4plot[,2], levels = c("Clinical Progressive Disease", "Stable Disease", "Partial Response", "Complete Response"))

pdf("4D.Paclitaxel.response.pdf", height=4, width=2.5)

p3 = ggplot(dat4plot, aes(x=response, y=Paclitaxel, fill=response)) +  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(shape=21, position = position_jitter(width = 0.3), size=0.5) + guides(colour=FALSE) +
  theme(axis.title.x=element_blank(), legend.title = element_blank() ) +
  labs(title=paste("CCLE-based Model\np = ", format(p,digits=3), sep=""), x="", y = "Predicted Response to Paclitaxel, TCGA-BRCA") +
  theme(legend.position = "none", plot.title = element_text(hjust=0.5), axis.text.x = element_text(angle = 45, vjust=0.7) ) +
  stat_summary(fun.data = give.n, geom = "text") +
  scale_x_discrete(labels=c("Clinical Progressive Disease" = "Clinical\nProgressive\nDisease", "Stable Disease" = "Stable\nDisease", "Partial Response" = "Partial\nResponse", "Complete Response"="Complete\nResponse"))
  
print(p3)
dev.off()
