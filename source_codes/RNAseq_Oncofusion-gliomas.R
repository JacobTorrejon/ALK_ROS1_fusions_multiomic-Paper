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
library(DESeq2)
library(DT)
library(factoextra)
library(ComplexHeatmap)
library(dendextend)
library(circlize)

##############################################################################################################
# 1 OPEN RAW COUNTS MATRICES AND SAMPLE PLAN
##############################################################################################################
### Raw counts data
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper/RNAseq/github/v2")
matrix1 <- read.delim("ALK_RNA_Feature_count.txt", header = TRUE)
# saveRDS(matrix1,"ALK_RNA_Feature_count.rds")
matrix2 <- read.delim("ROS1_RNA_feature_count.txt")
# saveRDS(matrix2,"ROS1_RNA_feature_count.rds")
RNA_rawcounts <- cbind(matrix1[,-c(2:6)], matrix2[,-c(1:6)])


### Load Sample Plan
Sample_Plan <- read.delim("Sample_Plan_RNAseq.txt")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
t <- c("Geneid",Sample_Plan$Replica_name)
colnames(RNA_rawcounts) <- t

##############################################################################################################
# 2 FILTERING TRANSCRIPTS
##############################################################################################################

#Remove transcript with zero in all samples
RNA_rawcounts_Nozero <- RNA_rawcounts[apply(RNA_rawcounts[,-c(1)],1,function(x) !all(x==0)),]

library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
Anno_Gene <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", 
                                "start_position", "end_position","band","transcript_biotype"),
                 mart=mart)

Anno_protcoding <- subset(Anno_Gene, Anno_Gene[,"transcript_biotype"] == "protein_coding")
list_protcoding <- Anno_protcoding$hgnc_symbol
RNA_rawcounts_filtered <- RNA_rawcounts_Nozero[RNA_rawcounts_Nozero$Geneid %in% list_protcoding , ]
row.names(RNA_rawcounts_filtered) <- RNA_rawcounts_filtered[,"Geneid"]

### Remove duplicates by taking the one with highest value
RNA_rawcounts_filtered$Mean <- rowMeans(RNA_rawcounts_filtered[,c(-1)])
RNA_rawcounts_filtered2 <- RNA_rawcounts_filtered[order(RNA_rawcounts_filtered$Mean, decreasing = TRUE), ]
RNA_rawcounts_filtered_Nodup <- RNA_rawcounts_filtered2[!duplicated(RNA_rawcounts_filtered2$Geneid), ]
row.names(RNA_rawcounts_filtered_Nodup) <- RNA_rawcounts_filtered_Nodup$Geneid
RNA_rawcounts_filtered_Nodup <- RNA_rawcounts_filtered_Nodup[,-c(1,ncol(RNA_rawcounts_filtered_Nodup))]
# write.table(RNA_rawcounts_filtered_Nodup,"Matrix_countsraw_Protein coding genes.txt",sep="\t",row.names=TRUE)

boxplot(RNA_rawcounts_filtered_Nodup[,])

##############################################################################################################
# 3 NORMALIZATION DESEQ2
##############################################################################################################

dds <- DESeqDataSetFromMatrix(countData=RNA_rawcounts_filtered_Nodup, DataFrame(condition=Sample_Plan$Oncofusion), ~ condition)
dds
dds <- estimateSizeFactors(dds)
dds <- dds[ rowSums(counts(dds)) > 0, ]

### rlog normalization
rld <- rlog(dds, blind=TRUE)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3, main="log2 normalized counts")
plot(assay(rld)[,1:2],
     pch=16, cex=0.3, main="rlog normalized counts")


plotPCA(rld, intgroup = c("condition"), ntop=5000)


d.rlog <- assay(rld)
boxplot(d.rlog, ylim=c(-2,20),main="rlog normalization")
gvar <- apply(d.rlog, 1, var)
mostvargenes <- order(gvar, decreasing=TRUE)[1:5000]
res_pca <- PCA(t(d.rlog[mostvargenes,]), ncp=3, graph=TRUE)

fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_ind(res_pca, axes=c(1,2), label="none",
             palette = c("ALK" = "#147082", "ROS1" = "#8F0A2E"),
             geom.ind = "point",
             pointshape = 21,
             pointsize = 2.5,
             fill.ind = as.factor(Sample_Plan$Oncofusion),
             title = "PCA (rlog Norm) - 5000 most variable genes")

fviz_pca_ind(res_pca, axes=c(1,2), label="none",
             palette = c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082",
                         "GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                         "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"),
             geom.ind = "point",
             pointshape = 21,
             pointsize = 2.5,
             fill.ind = as.factor(Sample_Plan$Sample_name),
             title = "PCA (rlog Norm) - 5000 most variable genes")



par(mar = c(13, 4, 4, 2) + 0.1);
boxplot(d.rlog,ylab="RNAseq level",las=2,main="rlog normalized counts")
mat <- as.matrix(d.rlog)



##########################################################
# 4 EXPLORATORY ANALYSIS (PCA)                             #
##########################################################

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
  
  
  xlab("PC1 (35.6 %)") +
  ylab("PC2 (14.9 %)") +
  ggtitle("PCA_RNAseq")



g1 <- g + scale_fill_manual(values = c("ALK" = "#147082", "ROS1" = "#8F0A2E"))

g1



# library(FactoMineR)
# library(factoextra)
# gvar <- apply(mat, 1, var)
# mostvargenes <- order(gvar, decreasing=TRUE)[1:nrow(mat)]
# res_pca <- PCA(t(mat[mostvargenes,]), ncp=10, graph=TRUE)
# fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main=("PCA_RNAseq"))





##########################################################
# 5. EXPLORATORY ANALYSIS (HC)                          #
##########################################################

### Most variable genes by condition / group
mat <- as.matrix(d.rlog)
Pval <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Sample_name)
  return(summary(trans)[[1]][1,5])
}

t.mat <- apply(mat, 1, Pval)
t.mat <- t.mat[t.mat <= 1e-20]
mat <- mat[names(t.mat),]



#ComplexHeatmap
SamplePlan2 <- as.data.frame(t(Sample_Plan))
Group <- SamplePlan2["Oncofusion",]
mat_G1 <- as.matrix(Group)
Cell_Line <- SamplePlan2["Sample_name",]
mat_CL <- as.matrix(Cell_Line)



mat <- t(scale(t(mat)))

colors <- c("black","black","black")
Nclust <- 5
column_tree = hclust(as.dist(1-cor(mat, method="pearson")), method = "ward.D")
column_tree = color_branches(column_tree, k = Nclust, col = colors)
# No imputation
# column_tree = hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
# column_tree = color_branches(column_tree, k = Nclust, col = colors)


#Cluster for rows
# column_tree2 = hclust(as.dist(1-cor(t(mat), method="pearson")), method = "ward.D")
column_tree2 = hclust(dist(mat, method = "euclidean"), method = "complete")
column_tree2 = color_branches(column_tree2, k = Nclust, col = colors)

ht <- Heatmap(mat, na_col = "white", 
              # col = colorRamp2(c(3, 5, 7, 9, 11), c("darkblue", "cyan", "white", "tomato", "darkred")),
              column_dend_reorder = FALSE,
              column_title = "HC TOP955 most variable genes",
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


### Subgroups extraction
# hc <- row_order(ht)
# 
# for (i in 1:length(hc)){
#   if (i == 1) {
#     clu <- t(t(row.names(mat[hc[[i]],])))
#     proteomic <- cbind(clu, paste("cluster_RNA_", i, sep=""))
#     colnames(proteomic) <- c("New_name", "Proteomic")
#   } else {
#     clu <- t(t(row.names(mat[hc[[i]],])))
#     clu <- cbind(clu, paste("cluster_RNA_", i, sep=""))
#     proteomic <- rbind(proteomic, clu)
#   }
# }
# 
# 
# proteomic
# colnames(proteomic) <- c("Gene","Cluster_RNA") 
# setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper")
# # write.table(proteomic,"list_Oncofusion_cluster_RNA_TOP955.txt",sep="\t",row.names=FALSE)



