##############################################################################################
### 0 LOAD LIBRARIES
##############################################################################################

library(stats)
library(plyr)
library(dplyr)
library(ggplot2)
library(tidyverse) 
library(tidyr)
library(ConsensusClusterPlus)
library(FactoMineR)
library(missMDA)
library(NMF)
library(ComplexHeatmap)
library(dendextend)
library(circlize)




##############################################################################################################
# 1 OPEN PHOSPHO MATRIx AND SAMPLE PLAN
##############################################################################################################

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper/Phospho/github/v2")
phospho <- readRDS("Phospho_matrix_ALK_ROS1_normalized.rds")


Sample_Plan <- read.csv("Sample_Plan_Phospho.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
Sample_ID <- Sample_Plan$Sample_ID

phospho <- phospho[,Sample_ID]
colnames(phospho) <- row.names(Sample_Plan)

par(mar = c(12, 4, 4, 2) + 0.1);
boxplot(phospho,las=2, ylab = "PTM level")

dev.off()

############################################################################################################################
# 2 REMOVE PTMs WITH MISSING VALUES
############################################################################################################################

nzero.phospho <- apply(phospho, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )

# Plot N.missing values vs N.quantified proteins
df.nzero.phospho <- data.frame(table(nzero.phospho))
x.df.nzero.phospho <- as.numeric(as.character(df.nzero.phospho[,1]))
x.df.nzero.phospho2 <- x.df.nzero.phospho*100/ncol(phospho)

y.df.nzero.phospho <- cumsum(df.nzero.phospho[,2])
plot(x.df.nzero.phospho2, y.df.nzero.phospho, pch = 16,  main = "", xlab = "% Missing values", 
     ylab = "Quantified proteins", cex=2, ylim = c(0,11000) )
grid()

## Remove PTMs with NAs
na.th = ncol(phospho)*0.01
phospho.10na <- phospho[which(nzero.phospho<=na.th),] 
boxplot(phospho.10na)
mean.centered.phospho.10na <- scale(phospho.10na, center = TRUE, scale = TRUE)

boxplot(mean.centered.phospho.10na)



##########################################################
# 3. EXPLORATORY ANALYSIS (PCA)                          #
##########################################################

# mat <- phospho.10na
mat <- mean.centered.phospho.10na

SamplePlan2 <- as.data.frame(t(Sample_Plan))
Group <- SamplePlan2["Oncofusion",]
mat_G1 <- as.matrix(Group)
Cell_Line <- SamplePlan2["Sample_name",]
mat_CL <- as.matrix(Cell_Line)

mat2 <- as.data.frame(rbind(mat_G1,mat_CL,mat))
mat3 <- as.data.frame(t(mat2))

mat_PCA <- mat3[,3:ncol(mat3)]
str(mat_PCA)
sapply(mat_PCA, class)
mat_PCA2 <- mat_PCA
mat_PCA2[] <- lapply(mat_PCA2, function(x) as.numeric(as.character(x)))
sapply(mat_PCA2, class)

myPr <- prcomp(mat_PCA2[,])
myPr
summary(myPr)
plot(myPr, type = "l")

biplot(myPr)
biplot(myPr, scale = 0)
str(myPr)
myPr$x

mat_PCA3 <- mat_PCA2
mat_PCA3 <- cbind(mat_PCA2,mat3$Oncofusion)
# mat_PCA3 <- cbind(mat_PCA2,mat3$Sample_name)

my_data_PCA <- cbind(mat_PCA3, myPr$x)
my_data_PCA2 <- as.data.frame(my_data_PCA)
names(my_data_PCA2)[names(my_data_PCA2) == 'mat3$Oncofusion'] <- 'Oncofusion'
# names(my_data_PCA2)[names(my_data_PCA2) == 'mat3$Sample_name'] <- 'Sample'


g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2))
                              , col = Oncofusion, fill = Oncofusion
)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(size = 5, shape = 21, col = "black") +
  theme_light(base_size=25) +
  
  
  xlab("PC1 (74.9 %)") +
  ylab("PC2 (9.2 %)") +
  ggtitle("PCA_Phospho")



g1 <- g + scale_fill_manual(values = c("ALK" = "#147082", "ROS1" = "#8F0A2E"))

g1

# library(FactoMineR)
# library(factoextra)
# gvar <- apply(mat, 1, var)
# mostvargenes <- order(gvar, decreasing=TRUE)[1:nrow(mat)]
# res_pca <- PCA(t(mat[mostvargenes,]), ncp=10, graph=TRUE)
# fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main=("PCA_Phospho"))
# 


##########################################################
# 4 EXPLORATORY ANALYSIS (HC)                        #
##########################################################

### Most variable proteins by condition / group
# mat <- phospho.10na
mat <- mean.centered.phospho.10na

Pval <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Sample_name)
  return(summary(trans)[[1]][1,5])
}

t.mat <- apply(mat, 1, Pval)
t.mat <- t.mat[t.mat <= 1e-25]
mat <- mat[names(t.mat),]



#ComplexHeatmap
SamplePlan2 <- as.data.frame(t(Sample_Plan))
Group <- SamplePlan2["Oncofusion",]
mat_G1 <- as.matrix(Group)
Cell_Line <- SamplePlan2["Sample_name",]
mat_CL <- as.matrix(Cell_Line)



mat <- t(scale(t(mat)))

colors <- c("black","black","black")
Nclust <- 2

#Cluster for columns
column_tree = hclust(as.dist(1-cor(mat, method="pearson")), method = "ward.D")
column_tree = color_branches(column_tree, k = Nclust, col = colors)
# column_tree = hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
# column_tree = color_branches(column_tree, k = Nclust, col = colors)


#Cluster for rows
# column_tree2 = hclust(as.dist(1-cor(t(mat), method="pearson")), method = "ward.D")
column_tree2 = hclust(dist(mat, method = "euclidean"), method = "complete")
column_tree2 = color_branches(column_tree2, k = Nclust, col = colors)

ht <- Heatmap(mat, na_col = "white", 
              # col = colorRamp2(c(3, 5, 7, 9, 11), c("darkblue", "cyan", "white", "tomato", "darkred")),
              column_dend_reorder = FALSE,
              column_title = "HC TOP1019 most variable PTMs",
              name = "Z score", 
              cluster_rows = column_tree2,
              cluster_columns = column_tree,
              column_split = 2,
              row_split = Nclust,
              show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 8), 
              column_names_gp = gpar(fontsize = 10), 
              top_annotation = HeatmapAnnotation(
                Oncofusion = mat_G1[1,],
                Cell_Line = mat_CL [1,],
                col = list(
                  Oncofusion = c("ALK" = "#147082", "ROS1" = "#8F0A2E"),
                  Cell_Line = c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082",
                                "GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                                "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")
                )
              ))

ht = draw(ht)



dev.off(which = dev.cur())


# ### Subgroups extraction
# hc <- row_order(ht)
# 
# for (i in 1:length(hc)){
#   if (i == 1) {
#     clu <- t(t(row.names(mat[hc[[i]],])))
#     proteomic <- cbind(clu, paste("cluster_PHOSPHO_", i, sep=""))
#     colnames(proteomic) <- c("New_name", "Proteomic")
#   } else {
#     clu <- t(t(row.names(mat[hc[[i]],])))
#     clu <- cbind(clu, paste("cluster_PHOSPHO_", i, sep=""))
#     proteomic <- rbind(proteomic, clu)
#   }
# }
# 
# 
# proteomic
# colnames(proteomic) <- c("Gene","Cluster_PHOSPHO") 
# setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper")
# # write.table(proteomic,"list_Oncofusion_cluster_PHOSPHO_TOP1019.txt",sep="\t",row.names=FALSE)




