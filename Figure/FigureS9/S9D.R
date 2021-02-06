setwd("/path/to/VAEN/Figure/FigureS9")
library("reshape2")
library("ggplot2")

load("../Figure3/3EF.data.RData")
rownames(TCGA.drug.by.cancer.mat) = rownames(obsd.drug.by.tissue.mat)

##################################################################
two.drugs.match = read.delim("../../DATA/drugs.match-2.txt", as.is=T, header=F)
two.drugs.match = two.drugs.match[order(two.drugs.match[,3]), ]

match(two.drugs.match[,3], rownames(TCGA.drug.by.cancer.mat)) -> shared.ii
shared.ii = shared.ii[length(shared.ii):1]

which( apply(obsd.drug.by.tissue.mat, 1, function(u)sum(abs(u) < 0.05)) > 8) -> ii
unique.ii = setdiff(ii, shared.ii)
unique.ii = unique.ii[length(unique.ii):1]
draw.ii = c( unique.ii, shared.ii )

############

draw.mat = TCGA.drug.by.cancer.mat[draw.ii,]
rownames(draw.mat) -> label
label[which(label == "Crizotinib")] = "PF2341066 | Crizotinib"
label[which(label == "NVP-TAE684")] = "TAE684 | NVP-TAE684"
label[which(label == "Palbociclib")] = "PD-0332991 | Palbociclib"
label[which(label == "Saracatinib")] = "AZD0530 | Saracatinib"
label[which(label == "Tanespimycin")] = "17-AAG | Tanespimycin"

rownames(draw.mat) = label
melt(draw.mat) -> mm
x = mm[,3]

###
which( abs(mm[,3]) < 1e-100) -> ii
mm[ii,3] = sign(x)[ii] * 1e-100
x = mm[, 3]
###

mm[,3] = -log10(abs(x)) * sign(x)

ggheatmap.3 <- ggplot(data = mm, aes(Var2, Var1, fill = value))+
 geom_tile(color = "grey")+
 scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 0, space = "Lab", name="t-test") +
 theme_minimal()+ 
 theme(axis.text.x.top = element_text(angle = 90, vjust = 0.5, size = 8, hjust = 0), axis.text.y = element_text(size = 8), axis.title.x=element_blank(), axis.title.y=element_blank() )+ scale_x_discrete(position="top")+
 coord_fixed()

pdf("FigureS9D.pdf", width=6.2, height=4.1)
print(ggheatmap.3)
dev.off()

