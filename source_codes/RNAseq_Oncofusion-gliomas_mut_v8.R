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

## We load the required packages
library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(OmnipathR)
library(dorothea)
library(tidyverse)
library(magrittr)


#OPEN DATA FILE----------------------------------------------------------------------------------------------
### 1. Raw couunts data
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/RNAseq")
matrix1 <- read.delim("ALK_RNA_Feature_count.txt", header = TRUE)
matrix2 <- read.delim("ROS1_feature_count_RNA.txt")
RNA_rawcounts <- cbind(matrix1[,-c(2:6)], matrix2[,-c(1:6)])



### 2. Load Sample Plan--------------------------------------------------------------------------------------------
Sample_Plan <- read.delim("Sample_Plan2.txt")
# Sample_Plan <- read.csv("Sample_Plan2.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
t <- c("Geneid",Sample_Plan$Replica_name)
colnames(RNA_rawcounts) <- t


### 3 Filter transcript-------------------------------------------------------------------------------------

#Remove transcript with zero in all samples
RNA_rawcounts_Nozero <- RNA_rawcounts[apply(RNA_rawcounts[,-c(1)],1,function(x) !all(x==0)),]

library(biomaRt)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
Anno_Gene <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "chromosome_name", 
                                "start_position", "end_position","band","transcript_biotype"),
                 # filters = c("chromosome_name", "start", "end"), 
                 # values=list(8, 1, 96000000),
                 mart=mart)

Anno_protcoding <- subset(Anno_Gene, Anno_Gene[,"transcript_biotype"] == "protein_coding")
list_protcoding <- Anno_protcoding$hgnc_symbol
RNA_rawcounts_filtered <- RNA_rawcounts_Nozero[RNA_rawcounts_Nozero$Geneid %in% list_protcoding , ]
row.names(RNA_rawcounts_filtered) <- RNA_rawcounts_filtered[,"Geneid"]

###Remove duplicates by taking the one with highest value
RNA_rawcounts_filtered$Mean <- rowMeans(RNA_rawcounts_filtered[,c(-1)])
RNA_rawcounts_filtered2 <- RNA_rawcounts_filtered[order(RNA_rawcounts_filtered$Mean, decreasing = TRUE), ]
RNA_rawcounts_filtered_Nodup <- RNA_rawcounts_filtered2[!duplicated(RNA_rawcounts_filtered2$Geneid), ]
row.names(RNA_rawcounts_filtered_Nodup) <- RNA_rawcounts_filtered_Nodup$Geneid
RNA_rawcounts_filtered_Nodup <- RNA_rawcounts_filtered_Nodup[,-c(1,ncol(RNA_rawcounts_filtered_Nodup))]
# write.table(RNA_rawcounts_filtered_Nodup,"Matrix_countsraw_Protein coding genes.txt",sep="\t",row.names=TRUE)

boxplot(RNA_rawcounts_filtered_Nodup[,])




### 4. Selecting one group-----------------------------------------------------------------------

# Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Group"] == "Oligo")
# # Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Group"] == "Oligo" | Sample_Plan[,"Group"] == "Normal")
# 
# list_sub <- row.names(Sample_Plan_sub)
# RNA_rawcounts_filtered_Nodup_sub <- RNA_rawcounts_filtered_Nodup[, colnames(RNA_rawcounts_filtered_Nodup) %in% list_sub]
# RNA_rawcounts_filtered_Nodup <- RNA_rawcounts_filtered_Nodup_sub


### 5. exploratory analysis-----------------------------------
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
# write.table(d.rlog,"RNAseq matrix_FSTL5_Normalized_filtered_protein coding genes.txt",sep="\t",row.names=TRUE)
boxplot(d.rlog, ylim=c(-2,20),main="rlog normalization")
gvar <- apply(d.rlog, 1, var)
mostvargenes <- order(gvar, decreasing=TRUE)[1:5000]
res_pca <- PCA(t(d.rlog[mostvargenes,]), ncp=3, graph=TRUE)

fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100))

fviz_pca_ind(res_pca, axes=c(1,2), label="none",
             palette = c("blue","red","green"),
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




### vst normalization
# vst <- vst(dds, blind=TRUE)
# plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
#      pch=16, cex=0.3, main="log2 normalized counts")
# plot(assay(vst)[,1:2],
#      pch=16, cex=0.3, main="vst normalized counts")
# 
# 
# plotPCA(vst, intgroup = c("condition"), ntop=5000)
# 
# 
# d.vst <- assay(vst)
# boxplot(d.vst, ylim=c(-1,20))
# gvar <- apply(d.vst, 1, var)
# mostvargenes <- order(gvar, decreasing=TRUE)[1:5000]
# res_pca <- PCA(t(d.vst[mostvargenes,]), ncp=3, graph=TRUE)
# 
# fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50))
# 
# fviz_pca_ind(res_pca, axes=c(1,2), label="none",
#              palette = c("blue","red","green"),
#              geom.ind = "point",
#              pointshape = 21,
#              pointsize = 2.5,
#              fill.ind = as.factor(Sample_Plan$Group),
#              title = "PCA (vst Norm) - 5000 most variable genes")
# 
# 
# 
# boxplot(d.vst,ylab="RNAseq level", main = "vst normalization")
# mat <- as.matrix(d.vst)


##########################################################
# PCA with subgroups------------------------------------ #
##########################################################

SamplePlan2 <- as.data.frame(t(Sample_Plan))
# SamplePlan3 <- as.data.frame(t(Sample_Plan2_sub))
Group <- SamplePlan2["Oncofusion",]
mat_G1 <- as.matrix(Group)
Cell_Line <- SamplePlan2["Sample_name",]
mat_CL <- as.matrix(Cell_Line)

proteomic9 <- as.data.frame(rbind(mat_G1,mat_CL,mat))

proteomic10 <- as.data.frame(t(proteomic9))

colnames(proteomic10[,1:4])

mat_PCA <- proteomic10[,3:ncol(proteomic10)]
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
mat_PCA3 <- cbind(mat_PCA2,proteomic10$Oncofusion)
# mat_PCA3 <- cbind(mat_PCA2,proteomic10$Sample_name)



my_data_PCA <- cbind(mat_PCA3, myPr$x)

my_data_PCA2 <- as.data.frame(my_data_PCA)

names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Oncofusion'] <- 'Oncofusion'
# names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Sample_name'] <- 'Sample'



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

library(FactoMineR)
library(factoextra)
gvar <- apply(mat, 1, var)
mostvargenes <- order(gvar, decreasing=TRUE)[1:nrow(mat)]
res_pca <- PCA(t(mat[mostvargenes,]), ncp=10, graph=TRUE)
fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main=("PCA_RNAseq"))

nrow(mat)



##############################################################################################################
# 6. TF ACTIVITY FROM FILTERED MATRIX                                                                        #
##############################################################################################################

# Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Sample_name"] == "GOPC.ROS1" | Sample_Plan[,"Sample_name"] == "CLIP1.ROS1")
# Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Sample_name"] == "KIF21A.ROS1" | Sample_Plan[,"Sample_name"] == "CLIP1.ROS1")
# Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Sample_name"] == "GOPC.ROS1" | Sample_Plan[,"Sample_name"] == "KIF21A.ROS1")
# Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Sample_name"] == "CLIP1.ROS1" | Sample_Plan[,"Sample_name"] == "CLIP1.ROS1_KD")
# Sample_Plan_sub <- subset(Sample_Plan, Sample_Plan[,"Sample_name"] == "KIF21A.ROS1_KD" | Sample_Plan[,"Sample_name"] == "KIF21A.ROS1")
# list_sub <- row.names(Sample_Plan_sub)
# rna_p <- mat[, list_sub]

# genes <- Diff_rna_filter2$Prot
# rna_p <- mat[row.names(mat) %in% genes , list_sub]




Pval <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Sample_name)
  return(summary(trans)[[1]][1,5])
}

rna_p <- mat
t.rna <- apply(rna_p, 1, Pval)
t.rna <- t.rna[t.rna <= 0.01]
rna_p <- rna_p[names(t.rna),]


# Pval2 <- function(x)
# {
#   trans <- aov(x ~ Sample_Plan$Kinase_Dead)
#   return(summary(trans)[[1]][1,5])
# }
# 
# 
# rna_p <- mat
# t.rna <- apply(rna_p, 1, Pval2)
# t.rna <- t.rna[t.rna <= 0.05]
# rna_p <- rna_p[names(t.rna),]
# 



# Remove NAs and set row names
counts <- as.matrix(rna_p)
head(counts)


design <- Sample_Plan[,c("Replica_name","Sample_name","Kinase_Dead","Oncofusion")]
colnames(design) <- c("sample","Cell_Line","Kinase_Dead","Oncofusion")
design


# OmnipathR::collectri()
net <- get_collectri(organism='human', split_complexes=FALSE)
net

# Run ulm
sample_acts <- run_ulm(mat=counts, net=net, .source='source', .target='target',
                       .mor='mor',  minsize = 5)

sample_acts



# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  filter(statistic == 'ulm') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


# Scale per sample
# sample_acts_mat <- scale(sample_acts_mat)

# Get top tfs by Pvalue
sample_acts_mat2 <- t(sample_acts_mat)

t.sample_acts <- apply(sample_acts_mat2, 1, Pval)
t.sample_acts <- t.sample_acts[t.sample_acts <= 1e-23]
sample_acts_mat2 <- sample_acts_mat2[names(t.sample_acts),]



# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("#3A63BE", "white", "#A93426"))(palette_length)

my_breaks <- c(seq(-2, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 2, length.out=floor(palette_length/2)))

# design_p <- design %>% arrange(Kinase_Dead,Oncofusion,Cell_Line)
# t_sample <- design_p$sample
# sample_acts_mat2 <- sample_acts_mat2[,t_sample]
# 
# metadata = as.data.frame(design_p[,c(3,4,2)])
# colnames(metadata) <- c("Kinase_Dead","Oncofusion","Cell_Line")
# row.names(metadata) <- design_p$sample
# annoCol<-list(Cell_Line=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4",
#                           "GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
#                           "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"),
#               Kinase_Dead=c("CTRL"="peachpuff","KD"="magenta"),
#               Oncofusion=c("ALK"="blue","ROS1"="red"))

design_p <- design %>% arrange(Oncofusion,Kinase_Dead,Cell_Line,)
t_sample <- design_p$sample
sample_acts_mat2 <- sample_acts_mat2[,t_sample]

metadata = as.data.frame(design_p[,c(4,3,2)])
colnames(metadata) <- c("Oncofusion","Kinase_Dead","Cell_Line")
row.names(metadata) <- design_p$sample
annoCol<-list(Cell_Line=c("CCDC88A.ALK"="darkblue","PPP1CB.ALK"="green4","GOPC.ROS1"="darkred","CLIP1.ROS1"="yellow3","KIF21A.ROS1"="orange3",
                          "CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","GOPC.ROS1_KD"="tomato","CLIP1.ROS1_KD"="yellow1",
                          "KIF21A.ROS1_KD"="orange1"),
              Kinase_Dead=c("CTRL"="peachpuff","KD"="magenta"),
              Oncofusion=c("ALK"="blue","ROS1"="red"))


Nclust <- 3
# Plot
pheatmap(sample_acts_mat2, border_color = NA, color=my_color, breaks = my_breaks,
         cluster_cols=F, cluster_rows=T, 
         scale="row",
         clustering_distance_cols="euclidean", clustering_method="ward.D",
         cutree_cols = Nclust,
         cutree_rows = Nclust,
         annotation_col=data.frame(metadata),
         annotation_colors = annoCol,
         show_rownames= TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize = 8,
         main = "Transcription Factors_Most variables across Cell Lines(TOP43)") 





dev.off()


## Boxplot
input_data2 <- cbind(design,sample_acts_mat)
set.seed(123)
TF <- "TADA2A"
data_Protein <- as.data.frame(input_data2[,c("condition",TF)])
colnames(data_Protein) <- c("Proteome", "Protein")
data_ProteinB <- data_Protein


# data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
#                                  levels = c("WNT","SHH","G3","G4"),
#                                  labels = c("WNT","SHH","G3","G4"))

par(mar = c(11, 5, 4, 1) + 0.1)
par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ Proteome, data = data_ProteinB, main = TF, boxwex = 0.5, las =2,
        lwd = 2, xlab = '', ylab = 'TF activity', col = c("blue","red","yellow","green"),
        # ylim=c(4,9.5),
        outline=FALSE)

stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=25, pch = 16, col = c("darkblue","darkred","yellow3","green3"))



#Protein selection by highest SD-------------------------------------------------------------------------------------------
library(multiClust)
#Probe Ranking
mat_500 <- probe_ranking(input=NULL, probe_number=500,
                         probe_num_selection="Fixed_Probe_Num", data.exp=mat, method="SD_Rank")
library(genefilter)
#Select proteins with highest IQR
mat500B <- varFilter(mat, var.func=IQR, var.cutoff=0.7, filterByQuantile=TRUE)


#Select proteins with highest SD
vars <- apply(mat, 1, sd)
mat500B <- mat[vars > quantile(vars, 0.95), ] 


mat <- as.matrix(mat500B)


mat <- as.matrix(d.rlog)

### Most variable proteins by condition / group
Pval <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Sample_name)
  return(summary(trans)[[1]][1,5])
}

t.mat <- apply(mat, 1, Pval)
t.mat <- t.mat[t.mat <= 1e-20]
mat <- mat[names(t.mat),]



#ComplexHeatmap-----------
library(ComplexHeatmap)

library(dendextend)
library(circlize)

SamplePlan2 <- as.data.frame(t(Sample_Plan))
# SamplePlan3 <- as.data.frame(t(Sample_Plan2_sub))
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
hc <- row_order(ht)

for (i in 1:length(hc)){
  if (i == 1) {
    clu <- t(t(row.names(mat[hc[[i]],])))
    proteomic <- cbind(clu, paste("cluster_RNA_", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(row.names(mat[hc[[i]],])))
    clu <- cbind(clu, paste("cluster_RNA_", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}


proteomic
colnames(proteomic) <- c("Gene","Cluster_RNA") 
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper")
# write.table(proteomic,"list_Oncofusion_cluster_RNA_TOP955.txt",sep="\t",row.names=FALSE)



#Monti ConsensusCluster-------------------------------------------------------------------------------------------------


rcc = ConsensusClusterPlus(mat,maxK=12,reps=100,pItem=0.8,pFeature=1,title="example2",distance="pearson",
                           innerLinkage="ward.D",finalLinkage ="ward.D",clusterAlg="hc", seed = 123)



# resICL = calcICL(rcc,title="example")
# rcc3 = ConsensusClusterPlus(mat,maxK=8,reps=100,pItem=0.8,pFeature=1,title="example3",distance="euclidean",clusterAlg="km")

Nclust <- 5
cc <- as.data.frame(rcc[[Nclust]][["consensusClass"]])
colnames(cc) <- c("Consensus_cluster")
cc2 <- as.data.frame(t(cc))


# data_1pept7A  <- data_1pept7[c("New_name","Methylome","DKFZ"),]
# cc3 <- rbind(cc2,data_1pept7A[, colnames(data_1pept7A)])
# cc4 <- as.data.frame(t(cc3))
# cc5 <- cc4 %>% arrange(Consensus_cluster,Methylome)

mat_cons <- as.data.frame(rcc[[Nclust]][["consensusMatrix"]])
colnames(mat_cons) <- colnames(cc2)
row.names(mat_cons) <- colnames(cc2)
mat_cons2 <- as.matrix(mat_cons)

SamplePlan3 <- as.data.frame(t(Sample_Plan2))
# SamplePlan3 <- as.data.frame(t(Sample_Plan2_sub))
Group1 <- SamplePlan3["Group",]
mat_G1 <- as.matrix(Group1)
Group2 <- SamplePlan3["Group_Oligo_Affy",]
mat_G2 <- as.matrix(Group2)
Group3 <- SamplePlan3["Group_Oligo_RNAseq",]
mat_G3 <- as.matrix(Group3)
Graveendel <- SamplePlan3["Graveendel_5000",]
mat_Grave <- as.matrix(Graveendel)

Grade <- SamplePlan3["Grade_Oligo",]
mat_Grade <- as.matrix(Grade)
IDH_mut <- SamplePlan3["IDHmt",]
mat_IDHmut <- as.matrix(IDH_mut)
IDH_mut1 <- SamplePlan3["IDH1_mutation",]
mat_IDHmut1 <- as.matrix(IDH_mut1)

Proteome <- SamplePlan3["Proteome",]
mat_Prot <- as.matrix(Proteome)
RNAseq <- SamplePlan3["RNAseq",]
mat_RNAseq <- as.matrix(RNAseq)


# Cluster for columns
colors <- c("black","black","black","black","black")
column_tree = hclust(as.dist(1-cor(mat_cons2, method="pearson")), method = "ward.D")
column_tree = color_branches(column_tree, k = Nclust, col = colors)
# column_tree = hclust(dist(t(mat_cons2), method = "euclidean"), method = "ward.D")
# column_tree = color_branches(column_tree, k = Nclust, col = colors)


# Cluster for rows
column_tree2 = hclust(as.dist(1-cor(t(mat_cons2), method="pearson")), method = "ward.D")
# column_tree2 = hclust(dist(mat_cons2, method = "euclidean"), method = "ward.D")
column_tree2 = color_branches(column_tree2, k = Nclust)



ht <- Heatmap(mat_cons2, na_col = "white", 
              # col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("white", "lightblue", "blue1", "blue2", "blue3")),
              # col = colorRamp2(c(0.8, 0.85, 0.9, 0.95, 1), c("darkblue", "lightblue", "beige", "red1", "darkred")),
              col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue", "lightblue", "beige", "red1", "darkred")),
              column_dend_reorder = FALSE,
              column_title = "Monti_PROGLIO_RNAseq_Oligo Group_5000 most variable genes_IQR (hc_Pearson/Ward.D)",
              name = "Consensus score", 
              cluster_rows = column_tree2,
              cluster_columns = column_tree,
              column_split = Nclust,
              row_split = Nclust,
              show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 5), 
              column_names_gp = gpar(fontsize = 6), 
              top_annotation = HeatmapAnnotation(
                RNAseq = mat_RNAseq[1,],
                Proteome = mat_Prot[1,],
                Group_WHO = mat_G1[1,],
                Affy_Oligo = mat_G2[1,],
                RNAseq_Oligo = mat_G3[1,],
                Graveendel = mat_Grave[1,],
                Histologic_Grade = mat_Grade[1,],
                IDH1_IDH2_genome = mat_IDHmut[1,],
                IDH1_proteome = mat_IDHmut1[1,],
                col = list(
                  RNAseq = c("GBM_RNA" = "green", "Oligo_RNA1" = "orange3", "Oligo_RNA2" = "greenyellow", "Oligo_RNA3" = "yellow", "Astro_RNA" = "blue"),
                  Proteome = c("GBM_prot" = "green", "Normal_prot" = "orange3", "Oligo_prot1" = "orange3", "Oligo_prot2" = "greenyellow", "Oligo_prot3" = "yellow", "Astro_prot" = "blue", "LGG_prot6" = "cyan",
                               "O_prot1" = "yellow",  "O_prot2" = "orange1",  "O_prot3" = "red", "N/A"="grey","Outlier"="grey"),
                  Group_WHO = c("Astro" = "blue","GBM" = "green", "Oligo" = "yellow", "Normal"="orange3", "N/A"="grey"),
                  Affy_Oligo = c("O1" = "yellow1","O2" = "orange", "O3" = "red", "N/A"="grey"),
                  RNAseq_Oligo = c("O1" = "yellow1","O2" = "orange", "O3" = "red", "No_Consensus"="grey", "N/A"="grey","N"="orange3"),
                  Graveendel = c("IGS_0" = "orange3","IGS_9" = "yellow2", "IGS_17" = "blue", "IGS_18"="green"),
                  Histologic_Grade = c("2"="yellow1", "3"="yellow4","N/A"="grey"),
                  IDH1_IDH2_genome = c("Oui"="Black","Non"="peachpuff","N/A"="grey"),
                  IDH1_proteome = c("Oui"="Black","Non"="peachpuff","N/A"="grey"))
              )
              # left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = colors))) 
              
)


ht = draw(ht)

dev.off()

# Monte Carlo Clustering-------------------------------
library(M3C)
library(NMF) # loading for a heatmap plotting function
library(gplots) # loading this for nice colour scale
library(ggsci) # more cool colours

pca(mat,legendtextsize = 10,axistextsize = 10,dotsize=2)

res <- M3C(mat, maxK = 4, method = 1, removeplots = FALSE, iters=100, cores=4, seed = 456, clusteralg = "hc")
# res$scores
# res$plots
# M3C

Nclust <- 2
M3C <- res[["realdataresults"]][[Nclust]][["consensus_matrix"]]
data <- res$realdataresults[[Nclust]]$ordered_data # this is the data
annon <- res$realdataresults[[Nclust]]$ordered_annotation
# pca(data,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)
colnames(M3C) <- colnames(data)
mat_M3C <- as.matrix(M3C)


SamplePlan3 <- as.data.frame(t(Sample_Plan2))
# SamplePlan3 <- as.data.frame(t(Sample_Plan2_sub))
M3C <- as.data.frame(M3C)
SamplePlan3 <- SamplePlan3[names(M3C)]

Group1 <- SamplePlan3["Group",]
mat_G1 <- as.matrix(Group1)


Prot_M3C <- res[["realdataresults"]][[Nclust]][["ordered_annotation"]]
str(Prot_M3C)
Prot_M3C[] <- lapply(Prot_M3C, function(x) as.character(as.factor(x)))
sapply(Prot_M3C,class)
Prot_M3C$New_name <- row.names(Prot_M3C)


# Correspondance cluster number and proteome subtype
# Prot_M3C1 <- replace(Prot_M3C[,], Prot_M3C[,] == "1", "Oligo_RNA1")
# Prot_M3C2 <- replace(Prot_M3C1[,], Prot_M3C1[,] == "2", "Oligo_RNA3")
# Prot_M3C2A <- replace(Prot_M3C2[,], Prot_M3C2[,] == "3", "Oligo_RNA2")
# Prot_M3C2B <- replace(Prot_M3C2A[,], Prot_M3C2A[,] == "4", "GBM_RNA")
# Prot_M3C2C <- replace(Prot_M3C2B[,], Prot_M3C2B[,] == "5", "Astro_RNA")


# Prot_M3C2G <- subset(Prot_M3C2C, select = -c(New_name))
# # Prot_M3C2G <- subset(Prot_M3C2A, select = -c(New_name))
# colnames(Prot_M3C2G) <- c("RNAseq")
# mat_RNA <- as.matrix(t(Prot_M3C2G))


ht <- Heatmap(mat_M3C, na_col = "white", 
              col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue", "lightblue", "beige", "red1", "darkred")),
              column_dend_reorder = FALSE,
              column_title = "M3C_PROGLIO_Oligo Group_RNAseq_5000 most variable (IQR) genes",
              name = "Consensus score", 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              # column_split = Nclust,
              # show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 5), 
              column_names_gp = gpar(fontsize = 6), 
              top_annotation = HeatmapAnnotation(
                # Treatment = mat_G1[1,],
                Tumor = mat_G1[1,],
                col = list(
                  # Treatment = c("Control"="blue","LDE225"="red"))
                  Tumor = c("Primary"="blue","Relapse"="red"))
              ))

ht = draw(ht)


proteomic5C <- mat_G1

# proteomic5D <- as.data.frame(t(rbind(mat_RNA,mat_Prot,mat_G1,mat_G2,mat_Grade,mat_IDHmut,mat_IDHmut1)))
# # write.table(proteomic5D,"RNAseq Classification_order by M3C clusters.txt",sep="\t",row.names=TRUE)
# 
# proteomic5E <- proteomic5D %>% arrange(row.names(proteomic5D))
# # write.table(proteomic5E,"RNAseq Classification_order by Sample_name.txt.txt",sep="\t",row.names=TRUE)



#Proteome matrix with proteome classification-----------------------


# data_1pept2B <- data_1pept2
# sapply(data_1pept2B,class)
# data_1pept2B[] <- lapply(data_1pept2B, function(x) as.factor(as.numeric(x)))
# sapply(data_1pept2B,class)

proteomic6 <- rbind(proteomic5C,d.rlog[, colnames(proteomic5C)])
sapply(proteomic6,class)

# saveRDS(proteomic6, "Gliomes_121tumors_Spectronaut_MyProMS norm_LFQ_16062021_withSamplePlan.rds")



# test1 <- data_1pept2["ME3","LGLIO_3334"]
# test2 <- proteomic6[c("Proteome","ME3"),"LGLIO_3334"]
# head(test1)
# head(test2)
# head(proteomic5C[,"LGLIO_3334"])





# PCA with subgroups------------------------------------

proteomic9 <- as.data.frame(rbind(proteomic5C,mat[, colnames(proteomic5C)]))
proteomic10 <- as.data.frame(t(proteomic9))


mat_PCA <- proteomic10 %>% select(-contains("RNAseq"),-contains("Proteome"),-contains("Group"),-contains("Group_Oligo")) 
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

# mat_PCA3 <- cbind(mat_PCA2,proteomic10$RNAseq)
# mat_PCA3 <- cbind(mat_PCA2,proteomic10$Proteome)
# mat_PCA3 <- cbind(mat_PCA2,proteomic10$Group)
mat_PCA3 <- cbind(mat_PCA2,proteomic10$Group_Oligo)

my_data_PCA <- cbind(mat_PCA3, myPr$x)

my_data_PCA2 <- as.data.frame(my_data_PCA)

# names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$RNAseq'] <- 'RNAseq'
# names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Proteome'] <- 'Proteomic'
# names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Group'] <- 'Group'
names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Group_Oligo'] <- 'Group_Oligo'


# g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2)), col = RNAseq, fill = RNAseq)) +
#   stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
#   geom_point(size = 4, shape = 21, col = "black")
# 
# # g1 <- g + scale_fill_manual(values = c("blue", "green","orange3","greenyellow","yellow"))
# g1 <- g + scale_fill_manual(values = c("orange3","greenyellow","yellow"))
# 
# g1 + geom_text(col = "black", size = 3, check_overlap = TRUE, vjust = 0, nudge_y = 0.4)



# g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2)), col = Proteomic, fill = Proteomic)) +
#   stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
#   geom_point(size = 4, shape = 21, col = "black")
# 
# g1 <- g + scale_fill_manual(values = c("blue","green", "grey","orange3","greenyellow","yellow","grey"))
# 
# g1 + geom_text(col = "black", size = 3, check_overlap = TRUE, vjust = 0, nudge_y = 0.4)



# g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2)), col = Group, fill = Group)) +
#   # stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
#   geom_point(size = 4, shape = 21, col = "black")
# 
# g1 <- g + scale_fill_manual(values = c("blue", "green","grey","orange3","yellow"))
# 
# g1 + geom_text(col = "black", size = 3, check_overlap = TRUE, vjust = 0, nudge_y = 0.4)



g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2)), col = Group_Oligo, fill = Group_Oligo)) +
  # stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(size = 4, shape = 21, col = "black")

g1 <- g + scale_fill_manual(values = c("grey","yellow","orange","red"))

g1 + geom_text(col = "black", size = 3, check_overlap = TRUE, vjust = 0, nudge_y = 0.4)



#PCA in 3D (FORMAT1)--------------------------------------------------------------------------------------
library(rgl)
library(car)


PC1 <- my_data_PCA2[,"PC1"]
PC2 <- my_data_PCA2[,"PC2"]
PC3 <- my_data_PCA2[,"PC3"]
Cluster <- as.factor(my_data_PCA2[,"RNAseq"])

scatter3d(x = PC1, y = PC2, z = PC3, groups = Cluster,
          surface=FALSE, grid = FALSE, ellipsoid = TRUE,
          surface.col = c("blue", "green3","orange3","greenyellow","yellow"),
          axis.col = c("black", "black", "black"),sphere.size=1.15, radius=1)


rgl.snapshot(filename = "PCA3D_histoGroup.png")


#PCA in 3D (FORMAT2)--------------------------------------------------------------------------------------
library(plot3D)

colors <- c("green", "blue")
colors <- colors[as.numeric(my_data_PCA2$Proteomic)]
p3d <- plot3d(my_data_PCA2$PC1, my_data_PCA2$PC2, my_data_PCA2$PC3, col= colors, size=1, type='s' 
)
text3d(my_data_PCA2$PC1, my_data_PCA2$PC2, my_data_PCA2$PC3+1,
       texts=c(rownames(my_data_PCA2)),col= colors , cex= 0.7, pos=3)

rgl.snapshot(filename = "PCA3D_RNAseq.png")





#UMAP-----------------------------------------------------------------------------------

library(umap)
set.seed(123)
umap <- umap(mat_PCA2,n_neighbors = 15L)


df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 RNAseq = proteomic10$RNAseq
                 # Proteome = proteomic10$Proteomic,
                 # Group = proteomic10$Group,
                 # Group_Oligo = proteomic10$Group_Oligo
                 )


t <- ggplot(df, aes(x, y, colour = RNAseq)) +
  geom_point(size = 6,pch = 16) +
  theme_light(base_size=20) +


  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("UMAP_PROGLIO_RNAseq")

t1 <- t + scale_colour_manual(values = c("blue", "green","orange3","greenyellow","yellow"))

t1


# t <- ggplot(df, aes(x, y, colour = Proteome)) +
#   geom_point(size = 6,pch = 16) +
#   theme_light(base_size=20) +
# 
# 
#   xlab("UMAP1") +
#   ylab("UMAP2") +
#   ggtitle("UMAP_LGG_hc")
# 
# t1 <- t + scale_colour_manual(values = c("green", "orange3","greenyellow","yellow","blue"))
# 
# t1


# df$Group <- factor(df$Group, exclude = NULL,
#                        levels = c("Astro","Cerveau","GBM","Oligo" , NA),
#                        labels = c("Astro","Cerveau","GBM","Oligo" , "N/A"))
# 
# 
# 
t <- ggplot(df, aes(x, y, colour = Group)) +
  geom_point(size = 6,pch = 16) +
  theme_light(base_size=20) +
  
  
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("UMAP_LGG_hc")

t1 <- t + scale_colour_manual(values = c("blue", "green","yellow","grey"))

t1


# df$Group_Oligo <- factor(df$Group_Oligo, exclude = NULL,
#                    levels = c("NR","O1","O2","O3", NA),
#                    labels = c("NR","O1","O2","O3","N/A"))
# 
# 
# 
# t <- ggplot(df, aes(x, y, colour = Group_Oligo)) +
#   geom_point(size = 6,pch = 16) +
#   theme_light(base_size=20) +
#   
#   
#   xlab("UMAP1") +
#   ylab("UMAP2") +
#   ggtitle("UMAP_LGG_hc") 
# 
# t1 <- t + scale_colour_manual(values = c("orange","blue","green","yellow","grey"))
# 
# t1







### Differential analysis-------------------------------------

dds <- DESeqDataSetFromMatrix(countData = RNA_rawcounts_filtered_Nodup,
                              colData = Sample_Plan,
                              design = ~ Oncofusion)



# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]

dds_Diff_Analysis <- DESeq(dds)
resultsNames(dds_Diff_Analysis)

# Extract results
res_Diff_Analysis <- results(dds_Diff_Analysis, 
                             contrast = c("Oncofusion",
                                          "ALK",
                                          "ROS1"
                             )
)

res_Diff_Analysis 

## To check if there is a difference between these two groups: pvalue distribution
ggplot(data = as.data.frame(res_Diff_Analysis), mapping = aes(x=pvalue))+ geom_histogram()

# Histogram of padj
ggplot(data = as.data.frame(res_Diff_Analysis), 
       mapping = aes(x=padj)
) + 
  geom_histogram()

#Density plot
ggplot(data = as.data.frame(res_Diff_Analysis),
       mapping=aes(x=padj)
)+
  geom_density()

res_Diff_Analysis2 <- as.data.frame(res_Diff_Analysis)

#ENHANCED VOLCANO PLOTS FC----------------------------------------------------------------
library(EnhancedVolcano)


Prot <- row.names(res_Diff_Analysis)
LogFC <- res_Diff_Analysis$log2FoldChange
Pvalue <- res_Diff_Analysis$padj
logPvalue <- -log10(Pvalue)

df <- as.data.frame(cbind(Prot,LogFC,Pvalue,logPvalue))
row.names(df) <- df[,"Prot"]
sapply(df, class)
# write.table(df,"Oncofusion-gliomas_ALK vs ROS1.txt",sep="\t",row.names=FALSE)

df[,2:ncol(df)] <- lapply(df[,2:ncol(df)], function(x) as.numeric(as.character(x)))
sapply(df, class)

df <- df[!is.na(df$LogFC),]
df <- df[!is.na(df$Pvalue),]

Rank_df <- df[order(df$LogFC, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_df,"RNAseq_Differential Analysis_Oncofusion gliomas_ALK vs ROS1.txt",sep="\t",row.names=FALSE)



IA2T2 <- df[(df$LogFC > 12),]
IA2T1 <- df[(df$LogFC < -3.5),]
IA2T1 <- IA2T1[(IA2T1$Pvalue < 0.00001),]
IA2T2 <- IA2T2[(IA2T2$Pvalue < 0.00001),]
IA2T0 <- df[(df$Pvalue < 1e-400),]


EnhancedVolcano(df,
                lab = rownames(df),
                x = 'LogFC',
                y = 'Pvalue',
                ylim = c(0,310),
                xlim = c(-6,17),
                # selectLab = c(rownames(IA2T2)),
                selectLab = c(rownames(IA2T2),rownames(IA2T1),rownames(IA2T0)),
                title = 'ALK1 vs ROS1 (Oncofusion driven gliomas)',
                pCutoff = 0.05,
                FCcutoff = 1,
                # pCutoff = 1e-40,
                # FCcutoff = 0.2,
                pointSize = 3.0,
                labSize = 3,
                col=c('grey', 'grey', 'grey', 'red4'),
                colAlpha = 1,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                colConnectors = 'black',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                max.overlaps= Inf
)



dev.off()


### Boxplots------------------

# Sample_Plan_sub2 <- as.data.frame(t(Sample_Plan_sub))
Sample_Plan_sub2 <- as.data.frame(t(Sample_Plan))

### TPM count
# input_data <- log2(RNA_rawcounts_filtered_Nodup+1)
# input_data2 <- cbind(Sample_Plan2,t(input_data))
# input_data2[,3:ncol(input_data2)] <- lapply(input_data2[,3:ncol(input_data2)], function(x) as.numeric(as.character(x)))
# sapply(input_data2, class)


# Normalized by DESEQ2
input_data <- as.data.frame(rbind(Sample_Plan_sub2,d.rlog[, colnames(Sample_Plan_sub2)]))
# input_data <- as.data.frame(rbind(Sample_Plan_sub2,d.vst[, colnames(Sample_Plan_sub2)]))
input_data2 <- as.data.frame(t(input_data))
input_data2[,6:ncol(input_data2)] <- lapply(input_data2[,6:ncol(input_data2)], function(x) as.numeric(as.character(x)))
sapply(input_data2, class)

set.seed(123)
Protein <- "PDPN"
data_Protein <- as.data.frame(input_data2[,c("Oncofusion",Protein)])
colnames(data_Protein) <- c("RNAseq", "Protein")
data_ProteinB <- data_Protein

data_ProteinB$RNAseq <- factor(data_ProteinB$RNAseq, exclude = NULL,
                               levels = c("ALK","ROS1"),
                               labels = c("ALK","ROS1"))

par(cex.lab=1.4, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ RNAseq, data = data_ProteinB, main = Protein, boxwex = 0.5,
        lwd = 2, xlab = 'Oncofusion', ylab = 'RNA level', col = c("blue","red1"), outline=FALSE)

stripchart(Protein ~ RNAseq, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=6, pch = 20, col = c("darkblue","darkred"))



#FOR SEVERAL PROTEINS---------------------------

input_data3 <- input_data2[,c("Group", "NTN1","UNC5A","UNC5B","UNC5C","UNC5D","DCC","NEO1")]
# write.table(input_data3,"NETRIN_RNAseq_PrimaryvsRelapse.txt",sep="\t",row.names=TRUE)
# write.table(input_data3,"NETRIN_RNAseq_LDE225.txt",sep="\t",row.names=TRUE)

sapply(input_data3,class)


var_list = combn(colnames(input_data3)[2:8], 1, simplify=FALSE)

par(mfrow=c(3,3))
# par(mfrow=c(15,15))

for (i in var_list){
  data_Protein <- as.data.frame(input_data3[,c("Group",i)])
  Protein <- i
  colnames(data_Protein) <- c("RNAseq", "Protein")
  data_ProteinB <- data_Protein
  data_ProteinB$RNAseq <- factor(data_ProteinB$RNAseq, exclude = NULL,
                                 levels = c("Control","LDE225"),
                                 labels = c("Control","LDE225"))

  # data_ProteinB$RNAseq <- factor(data_ProteinB$RNAseq, exclude = NULL,
  #                                levels = c("Primary","Relapse"),
  #                                labels = c("Primary","Relapse"))
  par(cex.lab=1.4, cex.axis=1.4, cex.main = 2.0)
  boxplot(Protein ~ RNAseq, data = data_ProteinB, main = Protein, boxwex = 0.5,
          lwd = 2, xlab = 'Tumor', ylab = 'Log2(TPM+1)', col = c("blue1","red1"), outline=FALSE)
  
  stripchart(Protein ~ RNAseq, vertical = TRUE, data = data_ProteinB,
             method = "jitter", add = TRUE, lwd=6, pch = 20, col = c("blue4","red4"))
  
  
}

dev.off()
