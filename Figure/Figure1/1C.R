setwd("/path/to/VAEN/Figures/Figure1")
library("tsne")

###################################################################################################

latent = read.table("/path/to/VAEN/result/1.CCLE.latent.tsv", as.is=T, header=T)

tsne(latent) -> tpc

sapply(rownames(latent), function(x){
	strsplit(x, split="\\.")[[1]][1] -> u
	strsplit(u, split="_")[[1]] -> v; v= v[-1]; paste(v, collapse="_")
}) -> tt
tissues = sort(unique(tt))

#> table(tt)
#tt
#                 AUTONOMIC_GANGLIA                      BILIARY_TRACT                               BONE 
#                                24                                 11                                 33 
#                            BREAST             CENTRAL_NERVOUS_SYSTEM                             CERVIX 
#                                53                                 75                                  8 
#                       ENDOMETRIUM                         FIBROBLAST HAEMATOPOIETIC_AND_LYMPHOID_TISSUE 
#                                28                                 35                                187 
#                            KIDNEY                    LARGE_INTESTINE                              LIVER 
#                                34                                 57                                 25 
#                              LUNG                         OESOPHAGUS                              OVARY 
#                               198                                 31                                 56 
#                          PANCREAS                             PLEURA                           PROSTATE 
#                                49                                  9                                  8 
#                              SKIN                        SOFT_TISSUE                            STOMACH 
#                                55                                 54                                 37 
#                           THYROID          UPPER_AERODIGESTIVE_TRACT                      URINARY_TRACT 
#                                12                                 33                                 36 


##### create the color vector
library("RColorBrewer")
c(brewer.pal(n = 9, name = "Paired"), brewer.pal(n = 8, name = "BrBG"), brewer.pal(n = 9, name = "Set3"), brewer.pal(n = 8, name = "RdBu"), brewer.pal(n = 8, name = "RdGy"), brewer.pal(n = 8, name = "PiYG"), brewer.pal(n = 8, name = "PuBu")) -> cc
cc = cc[1:length(tissues)]
color = c()
plot(tpc[,1], tpc[,2], pch=20, cex=.8, xlab="tSNE_1", ylab="tSNE_2")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20 )
	color = rbind(color, cbind(rownames(latent)[ii], rep(cc[k], length(ii)))) 
}
legend("topright", fill=cc, legend=tissues)

##### plot

pdf("1C.CCLE-tSNE.pdf", width=7, height=7)
plot(tpc[,1], tpc[,2], pch=20, cex=1.5, xlab="tSNE_1", ylab="tSNE_2", col="white")
for(k in 1:length(tissues)){
	which(tt == tissues[k]) -> ii
	points(tpc[ii,1], tpc[ii, 2], col=cc[k], pch=20, cex=1.5 )
}

plot(1,1)
legend("topright", fill=cc, legend=tissues, cex=.9)
dev.off()

