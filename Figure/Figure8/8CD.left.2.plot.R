setwd("/path/to/VAEN/Figure/Figure8")

source("../../code/unfactor.R")
source("../../code/multiplot.R")

give.n <- function(x){
   return(c(y = max(x) , label = length(x)))
}

########################################################################################################

load("THCA.BRAF.plot.RData")
tmp = unfactor(tmp)

which(tmp[,"grp"]=="Braf") -> ii
tmp[ii, "grp"] = "Braf-like"

which(tmp[,"grp"]=="Ras") -> ii
tmp[ii, "grp"] = "Ras-like"

p1 = ggplot(tmp, aes(x=Mut, y=PLX4720, fill=Mut)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
  ggtitle( paste("PLX4720, THCA\nStratified by mutation", sep="") ) + xlab("") + ylab("ActArea") + 
  theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5)) +
  theme(legend.position = "none") +
  stat_summary(fun.data = give.n, geom = "text", size=3)

p2 = ggplot(tmp, aes(x=grp, y=PLX4720, fill=Mut)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
  ggtitle( paste("PLX4720, THCA\nStratified by mutation and subgroup", sep="") ) + xlab("") + ylab("ActArea") + guides(fill=guide_legend()) + 
  stat_summary(fun.data = give.n, geom = "text", size=3) + 
  theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))

p3 = ggplot(tmp, aes(x=grp, y=PLX4720, fill=grp)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
  ggtitle( paste("PLX4720, THCA\nStratified by subgroup", sep="") ) + xlab("") + ylab("ActArea") + 
  stat_summary(fun.data = give.n, geom = "text", size=3) + 
  theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5)) +
  theme(legend.position = "none")


layout <- matrix(c(1,2,2,3), nrow = 1, byrow = T)
pdf("THCA.pdf", width=6, height=3)
multiplot(plotlist = list(p3, p2, p1), layout = layout)
dev.off()


########################################################################################################

load("SKCM.BRAF.plot.RData")
tmp = unfactor(tmp)

p1 = ggplot(tmp, aes(x=Mut, y=PLX4720, fill=Mut)) + geom_boxplot(position=position_dodge(0.8), outlier.shape = NA) + geom_point(size=1, stroke=.3, alpha=0.2, position = position_jitterdodge(dodge.width = 0.8)) +
  ggtitle( paste("PLX4720, SKCM\nStratified by mutation", sep="") ) + xlab("") + ylab("ActArea") + guides(fill=guide_legend()) + 
  stat_summary(fun.data = give.n, geom = "text", size=3) + 
  theme(text = element_text(size=8), axis.text.x=element_text(color="black", size = 8),axis.text.y=element_text(color="black", size = 8), plot.margin=unit(c(1,1,1,1), "mm"), plot.title = element_text(size = 8, hjust = 0.5))

pdf("SKCM.pdf", width=2, height=3)
print(p1)
dev.off()
