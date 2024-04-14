library(ggplot2)
library(plotly)
library(limma)
library(pheatmap)

use_interactive_graphics(interactive=F)

#How to subset according to FA class
lipid_classes = rowData(d_invitro)$Class %in% c("PS", "SM", "TAG")
groups = d_invitro$treatment %in% c("mcm", "control")
d_subset = d_invitro[complete.cases(assay(d_invitro[,1:9])),groups ]

#Quality check samples
plot_samples(d_invitro, type = "tic", log = TRUE)
#PLot molecules
plot_molecules(d_invitro, "cv", measure = "Area")
#Lipid class abundance
plot_lipidclass(d_invitro, "boxplot")

#PCA for both lines, check cell line after treaetment
mvaresults = mva(d_invitro, measure="Area", method="PCA")
plot_mva(mvaresults, color_by="treatment", components = c(1,2))

#PCoA for 215
mvaresults = mva(d_invitro[,10:18], measure="Area", method="PCA")
plot_mva(mvaresults, color_by="treatment", components = c(1,2))

#PCA for 253
mvaresults = mva(d_invitro[complete.cases(assay(d_invitro)),1:9], measure="Area", method="PCA")
plot_mva(mvaresults, color_by="treatment", components = c(1,2))


mvaresults = mva(d_invitro, method = "OPLS-DA", group_col = "cell_line")
plot_mva(mvaresults, color_by="cell_line", components = c(1,2))
plot_mva_loadings(mvaresults, color_by="Class", top.n=15)
top<-top_lipids(mvaresults, top.n=50)
top_molecules<-top$Molecule

#write a text file with top_molecules
write.table(top_molecules, file="top_molecules_eto_control.txt", sep="\t", quote=F, row.names=F, col.names=F)

heatmap(assay(d_invitro))
heatmap(assay(d_invitro[row.names(assay(d_invitro)) %in% top_molecules,]))



