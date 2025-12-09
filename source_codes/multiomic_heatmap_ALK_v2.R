### 1. LOAD PACKAGES

library(MOFA2)
library(ggplot2)
library(readr)
library(tidyverse)
library(pheatmap)
library(gplots)
library(gridExtra)
library(knitr)
library(base)
library(shiny)
library(DT)
library(ComplexHeatmap)
library(dendextend)
library(circlize)

##############################################################################################################
# 1. IMPORT MATRICES AND SAMPLE PLAN                                                                         #
##############################################################################################################
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/RNAseq")
rna <- readRDS("RNAseq matrix_ALK_Normalized_filtered_protein coding genes.rds")
# row.names(rna)<-paste(row.names(rna),"RNA",sep="_")

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Proteome")
prot <- readRDS("Proteome matrix_ALK_Normalized.rds")
# row.names(prot)<-paste(row.names(prot),"PROT",sep="_")

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024/phosphoMS ROS1 and ALK SerThr and Tyr/ALK")
phospho <- readRDS("Phospho matrix_ALK_Normalized.rds")

# Long_phospho_name <- row.names(phospho)
Short_phospho_name <- sub("\\.","-Phospho:",row.names(phospho))
# Names <- as.data.frame(cbind(Long_phospho_name,Short_phospho_name))
row.names(phospho) <- Short_phospho_name



# Sample Plan
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Proteome")
Sample_Plan <- read.csv("Sample_Plan2.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
meta <- Sample_Plan[1:16,]
names(meta)[names(meta) == 'Replica_name'] <- 'sample'
names(meta)[names(meta) == 'Sample_name'] <- 'Cell_Line'
rownames(meta) <- meta$sample

meta$Cell_Line <- as.factor(meta$Cell_Line)
meta$Oncofusion <- as.factor(meta$Oncofusion)




##############################################################################################################
# 3. KEEP THE SAME PATIENTS FOR ALL MATRICES                                                                 #
##############################################################################################################

t <- Reduce(intersect, list(colnames(rna), colnames(prot), colnames(phospho)))

rna_p <- rna[, t]
prot_p <- prot[, t]
phospho_p <- phospho[, t]
meta_p <- meta[t,]


###################################################
# 4. HEATMAP RNA
###################################################

### RNAseq
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
MOFA_F1 <- read.delim("Rank_weights_RNA_Factor1.txt")
MOFA_F1$Gene <- row.names(MOFA_F1)
MOFA_F1B <- MOFA_F1 %>% mutate(Gene = (gsub("_RNA", "", Gene)))


n <- 250 #500
MOFA_F1C <- MOFA_F1B[c(1:n,(nrow(MOFA_F1B)-(n-1)):nrow(MOFA_F1B)),]


F1_rna <- MOFA_F1C$Gene


MOFA_F2 <- read.delim("Rank_weights_RNA_Factor2.txt")
MOFA_F2$Gene <- row.names(MOFA_F2)
MOFA_F2B <- MOFA_F2 %>% mutate(Gene = (gsub("_RNA", "", Gene)))

# MOFA_F2C <- MOFA_F2B[c(1:n,(nrow(MOFA_F2B)-(n-1)):nrow(MOFA_F2B)),]
MOFA_F2C <- MOFA_F2B[c((nrow(MOFA_F2B)-(n-1)):nrow(MOFA_F2B)),]
F2_rna <- MOFA_F2C$Gene

F12_rna <- unique(c(F1_rna,F2_rna))


meta_MOFA <- meta_p %>% arrange(Cell_Line)
meta_MOFA <- meta_MOFA %>%
  mutate(Split_Prot = case_when(Cell_Line == 'CCDC88A.ALK'  ~ 1,
                                Cell_Line == 'CCDC88A.ALK_KD'  ~ 2,
                                Cell_Line == 'PPP1CB.ALK'  ~ 3,
                                Cell_Line == 'PPP1CB.ALK_KD'  ~ 4
  ))   
t_MOFA <- row.names(meta_MOFA)
rna_MOFA <- rna[F12_rna,t_MOFA]

meta_MOFA2 <- as.data.frame(t(meta_MOFA))
mat_G1 <-  as.matrix(meta_MOFA2["Cell_Line",]) 

mat_rna <- t(scale(t(rna_MOFA)))
Split <- meta_MOFA$Split_Prot


### Clustering
colors <- c("black","black","black" ,"black")
Nclust <- 4 #5

# #Cluster for columns
# # column_tree = hclust(as.dist(1-cor(mat, method="pearson")), method = "ward.D")
# # column_tree = color_branches(column_tree, k = Nclust, col = colors)
# # column_tree = hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
# # column_tree = color_branches(column_tree, k = Nclust, col = colors)
# 
# 
#Cluster for rows
column_tree2 = hclust(dist(mat_rna, method = "euclidean"), method = "complete")
column_tree2 = color_branches(column_tree2, k = Nclust,col = colors)

ht1 <- Heatmap(mat_rna,
               col = colorRamp2(c(-1.5, 0, 1.5), c("#0028A5", "#C2C2C2", "#FF2600")),
               # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               # rect_gp = gpar(col = "white", lwd = 0.1),
               heatmap_legend_param = list(nrow = 1),
               row_title = "RNAseq",
               column_title = "Multiomic_Heatmap_ALK",
               name = "Z-score", 
               column_split = Split,
               # cluster_rows = FALSE,
               cluster_columns = FALSE,
               cluster_rows = column_tree2,
               row_split = Nclust,
               # cluster_columns = column_tree,
               show_row_names = FALSE, 
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 10), 
               column_names_gp = gpar(fontsize = 5), 
               top_annotation = HeatmapAnnotation(
                 Cell_Line = mat_G1[1,],
                 col = list(
                   Cell_Line = c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082")
                 )                             
               ))




ht = draw(ht1)



### Subgroups extraction
hc <- row_order(ht)

for (i in 1:length(hc)){
  if (i == 1) {
    clu <- t(t(row.names(mat_rna[hc[[i]],])))
    proteomic <- cbind(clu, paste("cluster_RNA_", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(row.names(mat_rna[hc[[i]],])))
    clu <- cbind(clu, paste("cluster_RNA_", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}

proteomic
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
# write.table(proteomic,"list_MOFA_cluster_RNA_k4.txt",sep="\t",row.names=FALSE)



###################################################
# 5. HEATMAP PROT
###################################################

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
MOFA_F1 <- read.delim("Rank_weights_Prot_Factor1.txt")
MOFA_F1$Gene <- row.names(MOFA_F1)
MOFA_F1B <- MOFA_F1 %>% mutate(Gene = (gsub("_Proteome", "", Gene)))


n <- 230
MOFA_F1C <- MOFA_F1B[c(1:n,(nrow(MOFA_F1B)-(n-1)):nrow(MOFA_F1B)),]


F1_rna <- MOFA_F1C$Gene


MOFA_F2 <- read.delim("Rank_weights_Prot_Factor2.txt")
MOFA_F2$Gene <- row.names(MOFA_F2)
MOFA_F2B <- MOFA_F2 %>% mutate(Gene = (gsub("_Proteome", "", Gene)))

# MOFA_F2C <- MOFA_F2B[c(1:n,(nrow(MOFA_F2B)-(n-1)):nrow(MOFA_F2B)),]
MOFA_F2C <- MOFA_F2B[c((nrow(MOFA_F2B)-(n-1)):nrow(MOFA_F2B)),]
F2_rna <- MOFA_F2C$Gene

F12_rna <- unique(c(F1_rna,F2_rna))


meta_MOFA <- meta_p %>% arrange(Cell_Line)
t_MOFA <- row.names(meta_MOFA)
prot_MOFA <- prot[F12_rna,t_MOFA]

meta_MOFA2 <- as.data.frame(t(meta_MOFA))
mat_G1 <-  as.matrix(meta_MOFA2["Cell_Line",]) 

mat_prot <- t(scale(t(prot_MOFA)))

### Clustering
colors <- c("black","black","black")
Nclust <- 4

# #Cluster for columns
# # column_tree = hclust(as.dist(1-cor(mat, method="pearson")), method = "ward.D")
# # column_tree = color_branches(column_tree, k = Nclust, col = colors)
# # column_tree = hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
# # column_tree = color_branches(column_tree, k = Nclust, col = colors)
# 
# 
#Cluster for rows
column_tree2 = hclust(dist(mat_prot, method = "euclidean"), method = "complete")
column_tree2 = color_branches(column_tree2, k = Nclust,col = colors)

ht2 <- Heatmap(mat_prot,
               col = colorRamp2(c(-1.5, 0, 1.5), c("#0028A5", "#C2C2C2", "#FF2600")),
               # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               # rect_gp = gpar(col = "white", lwd = 0.1),
               heatmap_legend_param = list(nrow = 1),
               row_title = "Proteome",
               column_title = NULL,
               name = "Z-score", 
               column_split = Split,
               # cluster_rows = FALSE,
               cluster_columns = FALSE,
               cluster_rows = column_tree2,
               row_split = Nclust,
               # cluster_columns = column_tree,
               show_row_names = FALSE, 
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 10), 
               column_names_gp = gpar(fontsize = 5), 
               # top_annotation = HeatmapAnnotation(
               #   Cell_Line = mat_G1[1,],
               #   col = list(
               #     Cell_Line = c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")
               #   )                             
               # )
               )



ht = draw(ht2)

### Subgroups extraction
hc <- row_order(ht)

for (i in 1:length(hc)){
  if (i == 1) {
    clu <- t(t(row.names(mat_prot[hc[[i]],])))
    proteomic <- cbind(clu, paste("cluster_PROT_", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(row.names(mat_prot[hc[[i]],])))
    clu <- cbind(clu, paste("cluster_PROT_", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}

proteomic
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
# write.table(proteomic,"list_MOFA_cluster_PROT_k4.txt",sep="\t",row.names=FALSE)

###################################################
# 6. HEATMAP PHOSPHO
###################################################
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
MOFA_F1 <- read.delim("Rank_weights_PTM_Factor1.txt")
MOFA_F1$Gene <- row.names(MOFA_F1)
MOFA_F1B <- MOFA_F1 %>% mutate(Gene = (gsub("_Phospho", "", Gene)))


n <- 500 #750
MOFA_F1C <- MOFA_F1B[c(1:n,(nrow(MOFA_F1B)-(n-1)):nrow(MOFA_F1B)),]


F1_rna <- MOFA_F1C$Gene


MOFA_F2 <- read.delim("Rank_weights_PTM_Factor2.txt")
MOFA_F2$Gene <- row.names(MOFA_F2)
MOFA_F2B <- MOFA_F2 %>% mutate(Gene = (gsub("_Phospho", "", Gene)))

# MOFA_F2C <- MOFA_F2B[c(1:n,(nrow(MOFA_F2B)-(n-1)):nrow(MOFA_F2B)),]
MOFA_F2C <- MOFA_F2B[c((nrow(MOFA_F2B)-(n-1)):nrow(MOFA_F2B)),]
F2_rna <- MOFA_F2C$Gene

F12_rna <- unique(c(F1_rna,F2_rna))



meta_MOFA <- meta_p %>% arrange(Cell_Line)
t_MOFA <- row.names(meta_MOFA)
phospho_MOFA <- phospho[F12_rna,t_MOFA]


meta_MOFA2 <- as.data.frame(t(meta_MOFA))
mat_G1 <-  as.matrix(meta_MOFA2["Cell_Line",]) 

mat_phospho <- t(scale(t(phospho_MOFA)))


### Clustering
colors <- c("black","black","black")
Nclust <- 4

# #Cluster for columns
# # column_tree = hclust(as.dist(1-cor(mat, method="pearson")), method = "ward.D")
# # column_tree = color_branches(column_tree, k = Nclust, col = colors)
# # column_tree = hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
# # column_tree = color_branches(column_tree, k = Nclust, col = colors)
# 
# 
#Cluster for rows
column_tree2 = hclust(dist(mat_phospho, method = "euclidean"), method = "complete")
column_tree2 = color_branches(column_tree2, k = Nclust,col = colors)

ht3 <- Heatmap(mat_phospho,
               col = colorRamp2(c(-1.5, 0, 1.5), c("#0028A5", "#C2C2C2", "#FF2600")),
               # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               # rect_gp = gpar(col = "white", lwd = 0.1),
               heatmap_legend_param = list(nrow = 1),
               row_title = "Phospho",
               column_title = NULL,
               name = "Z-score", 
               column_split = Split,
               # cluster_rows = FALSE,
               cluster_columns = FALSE,
               cluster_rows = column_tree2,
               row_split = Nclust,
               # cluster_columns = column_tree,
               show_row_names = FALSE, 
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 10), 
               column_names_gp = gpar(fontsize = 5), 
               # top_annotation = HeatmapAnnotation(
               #   Cell_Line = mat_G1[1,],
               #   col = list(
               #     Cell_Line = c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")
               #   )                             
               # )
               )

ht = draw(ht3)

hc <- row_order(ht)

for (i in 1:length(hc)){
  if (i == 1) {
    clu <- t(t(row.names(mat_phospho[hc[[i]],])))
    proteomic <- cbind(clu, paste("cluster_PHOSPHO_", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(row.names(mat_phospho[hc[[i]],])))
    clu <- cbind(clu, paste("cluster_PHOSPHO_", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}

proteomic

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
# write.table(proteomic,"list_MOFA_cluster_PHOSPHO_k4.txt",sep="\t",row.names=FALSE)


ht_list2 = ht1 %v% ht2 %v% ht3
draw(ht_list2)

# ht_list3 = ht1 + ht2 + ht3
# draw(ht_list3)

