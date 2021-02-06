setwd("/path/to/VAEN/Figure/FigureS9")
library("reshape2")
library("ggplot2")

load("../Figure3/3BCD.data.RData")

##################################################################
TCGA.drug.by.cancer.mat = TCGA.drug.by.cancer.mat[nrow(TCGA.drug.by.cancer.mat):1, ]
rownames(TCGA.drug.by.cancer.mat) -> ss
ss[which(ss=="X17.AAG")] = "17-AAG"
gsub("\\.", "-", ss) -> new.ss

new.ss[which(new.ss == "PD-0332991")] = "PD-0332991 | Palbociclib"

rownames(TCGA.drug.by.cancer.mat) = new.ss
draw.mat = TCGA.drug.by.cancer.mat

melt(draw.mat) -> mm
x = mm[,3]

which( abs(mm[,3]) < 1e-100) -> ii
mm[ii,3] = sign(x)[ii] * 1e-100
x = mm[, 3]

mm[,3] = -log10(abs(x)) * sign(x)

ggheatmap.3 <- ggplot(data = mm, aes(Var2, Var1, fill = value))+
 geom_tile(color = "grey")+
 scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 0, space = "Lab", name="t-test") +
 theme_minimal()+ 
 theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 0), axis.text.y = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank() )+ scale_x_discrete(position="top")+
 coord_fixed()

pdf("FigureS9C.pdf", width=6, height=4.5)
print(ggheatmap.3)
dev.off()

##################################################################

