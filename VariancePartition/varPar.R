library('variancePartition')
library(stringr)
library(ggplot2)
# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
#load("/VariancePartition/tutorial/varPartData.rda")

#form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)
#varPart <- fitExtractVarPartModel( geneExpr, form, info )
#vp <- sortCols( varPart )
#varPart@names
#plotPercentBars( vp[1:10,] )
#plotVarPart( vp )

##################################### START HERE #########################
labels <- read.csv("/VariancePartition/input_data/labels.csv")
labels$AD <- as.character(labels$AD)
labels$SCZ <- as.character(labels$SCZ)
#labels$SubID <- as.factor(labels$SubID)
labels$PLAQUE <- replace(labels$PLAQUE, is.na(labels$PLAQUE), 0)

ct_edges <- read.csv("/VariancePartition/input_data/ct-specific_tf-tg-links.csv")
ct_edges <- as.data.frame(ct_edges)
ct_edges_row <- ct_edges$X
ct_edges<-ct_edges[,-c(1)]
dim(ct_edges)
ct_edges <- as.matrix(ct_edges)
ct_edges_row <- str_wrap(ct_edges_row, width = 4)
row.names(ct_edges) <- ct_edges_row


form_1 <- ~ PLAQUE + Age +  (1|Ethinicity) + (1|Sex) + (1|AD)   #+ ~ PLAQUE

varPart_1 <- fitExtractVarPartModel( ct_edges, form_1, labels )
vp_1 <- sortCols( varPart_1 )

plotPercentBars( vp_1[1:20,] )
plotVarPart( vp_1 )

fig1 <- plotVarPart( vp_1 , label.angle = 45, main = "ct: TF-TG links", )
ggsave("/VariancePartition/plots/ct-specific_tf-tg-links_violin.png", fig1, dpi=300, width = 4, height = 3,)

fig2 <- plotPercentBars( vp_1[1:20,] )
ggsave("/VariancePartition/plots/ct-specific_tf-tg-links_VarPar_bar.png", fig2, dpi=300, width = 3, height = 3,)

write.csv(vp_1, file = "/VariancePartition/plots/ct-specific_tf-tg-links_VarPar.csv")
