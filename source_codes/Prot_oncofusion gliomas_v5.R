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

setwd("~/Documents/Curie_Jacob/Proteomic//Oncofusion driven gliomas/Proteome")

#OPEN DATA FILE----------------------------------------------------------------------------------------------
# Open matrix
matrix1 <- read.delim("ALK_TotalProteinTotalProteinQuantMatrix.txt", header = TRUE)
row.names(matrix1) <- matrix1$Accession
boxplot(matrix1[,5:10])

matrix2 <- read.delim("ROS1_PhosphoAnalysisTotalProteinQuantMatrix.txt")
row.names(matrix2) <- matrix2$Accession
boxplot(matrix2[,5:10])


matrix <- merge(x= matrix1[,-c(1)], y=matrix2[,-c(1)], by='row.names', all = TRUE)
list_prot <- matrix$Row.names

Annotation <- read.delim("p27310_Homosapiens_GRCh38_mapping_Ensp_Gene_transcript_symbol.txt")
Annotation$Proteinid <- substring(Annotation$Proteinid, 2)
Annotation$gene_symbol <- substring(Annotation$gene_symbol, 13)

Annotation_filtered <- Annotation[Annotation$Proteinid %in% list_prot , ]
row.names(Annotation_filtered) <- Annotation_filtered$Proteinid 
Annotation_filtered <- Annotation_filtered[list_prot,]

row.names(matrix) <- make.names(Annotation_filtered$gene_symbol, unique = TRUE)
boxplot(matrix[,2:35])

data_1pept <- matrix[,-c(1)]

data_1pept1 <- as.data.frame(scale(data_1pept, center = TRUE, scale = TRUE))
boxplot(data_1pept1)

data_1pept1["MYC",]


#Checking %NAs
sum(is.na(data_1pept))
sum(!is.na(data_1pept))
sum(is.na(data_1pept))*100/(ncol(data_1pept)*nrow(data_1pept))

# Sample_Plan <- read.delim("Sample_Plan.txt")
Sample_Plan <- read_csv("Sample_Plan2.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name


#DIAGNOSIS PLOTS FOR PCA IMPUTATION---------------------------------------------------------------------------------------------------

data_1pept2 <- data_1pept1
colnames(data_1pept2) <- row.names(Sample_Plan)
sapply(data_1pept2, class)

par(mar = c(12, 4, 4, 2) + 0.1);
boxplot(data_1pept2,las=2, ylab = "Protein level")

dev.off()

#Load data
prot <- as.matrix(data_1pept2)  ## Proteins in rows, samples in columns 
nzero.prot <- apply(prot, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )

# Plot N.missing values vs N.quantified proteins
df.nzero.prot <- data.frame(table(nzero.prot))
x.df.nzero.prot <- as.numeric(as.character(df.nzero.prot[,1]))
x.df.nzero.prot2 <- x.df.nzero.prot*100/ncol(prot)

y.df.nzero.prot <- cumsum(df.nzero.prot[,2])
plot(x.df.nzero.prot2, y.df.nzero.prot, pch = 16,  main = "", xlab = "% Missing values", 
     ylab = "Quantified proteins", cex=2, ylim = c(0,8500) )
grid()

## Create a matrix with at most 10 NAs per protein
na.th = ncol(prot)*0.1
prot.10na <- prot[which(nzero.prot<=na.th),] 

boxplot(prot.10na)

mean.centered.prot.10na <- scale(prot.10na, center = TRUE, scale = TRUE)

boxplot(mean.centered.prot.10na)
# res.comp = imputePCA(t(prot.10na),ncp=6,  threshold = 1e-06, maxiter = 1000)
# res.pca = PCA(res.comp$completeObs)
# prot.10na_imput <- res.comp$completeObs
# 
# 
# 
# #### Check PCA proteme data
# n.prot <- apply(prot.10na,2, function(x) sum(! is.na(x)))
# pca.prot <- PCA(scale(prot.10na_imput, scale=F ),scale.unit=TRUE,graph=F )
# 
# df.pca.prot <- data.frame(PC1 = pca.prot$ind$coord[, 1], PC2 = pca.prot$ind$coord[, 2], PC3 = pca.prot$ind$coord[, 3], Nprot = n.prot)
# 
# 
# p1 <- ggplot(data = df.pca.prot, aes(x = PC1, y = PC2, fill=Nprot)) +
#   geom_point(aes(size=Nprot), pch=21, colour="Black") +
#   theme_bw()
# ggExtra::ggMarginal(p1, type = "histogram", colour = "lightgray")
# 
# df.pca.prot %>%
#   gather(-Nprot, key = "var", value = "value") %>% 
#   ggplot(aes(x = value, y = Nprot, color=Nprot)) +
#   geom_point() +
#   ylim(0,7000)+
#   facet_wrap(~ var, scales = "free") +
#   theme_bw()
# 
# #### Check PCA proteme data after re-centering proteome matrix
# par(cex.axis=0.8, mar = c(12, 4, 4, 2) + 0.1);
# boxplot(prot.10na)  # Before re-centering
# mean.centered.prot.10na <- scale(prot.10na, center = TRUE, scale = FALSE)
# boxplot(mean.centered.prot.10na,las=2,ylab=c("Log10(LFQ)")) # After re-centering
# boxplot(prot.10na,las=2,ylab=c("Log10(LFQ)"))
# 
# res.comp.mean.centered = imputePCA(t(mean.centered.prot.10na),ncp=6,  threshold = 1e-06, maxiter = 1000)
# res.pca.mean.centered = PCA(res.comp.mean.centered$completeObs)
# mean.centered.prot.10na_imput <- res.comp.mean.centered$completeObs
# 
# 
# pca.prot.after <- PCA(scale(mean.centered.prot.10na_imput, scale=F ),scale.unit=TRUE,graph=F )
# 
# df.pca.prot.after <- data.frame(PC1 = pca.prot.after$ind$coord[, 1], PC2 = pca.prot.after$ind$coord[, 2], PC3 = pca.prot.after$ind$coord[, 3], PC4 = pca.prot.after$ind$coord[, 4], PC5 = pca.prot.after$ind$coord[, 5], Nprot = n.prot)
# 
# p2 <- ggplot(data = df.pca.prot.after, aes(x = PC1, y = PC2, fill=Nprot)) +
#   geom_point(aes(size=Nprot), pch=21, colour="Black") +
#   theme_bw()
# ggExtra::ggMarginal(p2, type = "histogram", colour = "lightgray")
# # p2
# 
# df.pca.prot.after %>%
#   gather(-Nprot, key = "var", value = "value") %>% 
#   ggplot(aes(x = value, y = Nprot, color=Nprot)) +
#   geom_point() +
#   ylim(0,4200)+
#   facet_wrap(~ var, scales = "free") +
#   theme_bw()
# 
# 
# 
# #Remove outliers
# mean.centered.prot.10na_imput_outliers <- data.frame(mean.centered.prot.10na_imput,Nprot = n.prot)
# mean.centered.prot.10na_imput_outliers2 <- mean.centered.prot.10na_imput_outliers[(mean.centered.prot.10na_imput_outliers$Nprot > 4730),]
# mean.centered.prot.10na_imput_outliers3 <- subset(mean.centered.prot.10na_imput_outliers2, select=-c(Nprot))
# 
# n.prot_out <- subset(mean.centered.prot.10na_imput_outliers2, select=c(Nprot))
# 
# 
# pca.prot.after_out <- PCA(scale(mean.centered.prot.10na_imput_outliers3, scale=F ),scale.unit=TRUE,graph=F )
# 
# df.pca.prot.after_out <- data.frame(PC1 = pca.prot.after_out$ind$coord[, 1], PC2 = pca.prot.after_out$ind$coord[, 2], PC3 = pca.prot.after_out$ind$coord[, 3], Nprot = n.prot_out)
# 
# p2 <- ggplot(data = df.pca.prot.after_out, aes(x = PC1, y = PC2, fill=Nprot)) +
#   geom_point(aes(size=Nprot), pch=21, colour="Black") +
#   theme_bw()+
#   guides(size = guide_legend(reverse = TRUE))
# ggExtra::ggMarginal(p2, type = "histogram", colour = "lightgray")
# # p2
# 
# df.pca.prot.after_out %>%
#   gather(-Nprot, key = "var", value = "value") %>% 
#   ggplot(aes(x = value, y = Nprot, color=Nprot)) +
#   geom_point() +
#   ylim(4000,5000)+
#   facet_wrap(~ var, scales = "free") +
#   theme_bw()



# mat <- t(prot.10na_imput)
# mat <- t(mean.centered.prot.10na_imput)
# mat <- t(mean.centered.prot.10na_imput_outliers3)

# write.table(mat,"Proteome_matrix_DDA_Osteosarcoma_LFQ_3pept_FilterNA_PCAimputation.txt",sep="\t",row.names=TRUE)

##########################################################
# PCA with subgroups------------------------------------ #
##########################################################

# mat <- prot.10na
mat <- mean.centered.prot.10na

SamplePlan2 <- as.data.frame(t(Sample_Plan))
# SamplePlan3 <- as.data.frame(t(Sample_Plan2_sub))
Group <- SamplePlan2["Oncofusion",]
mat_G1 <- as.matrix(Group)
Cell_Line <- SamplePlan2["Sample_name",]
mat_CL <- as.matrix(Cell_Line)

proteomic9 <- as.data.frame(rbind(mat_G1,mat_CL,mat))

proteomic10 <- as.data.frame(t(proteomic9))

colnames(proteomic10[,1:4])

mat_PCA <- proteomic10 %>% select(-contains("Oncofusion")) %>% select(-contains("Sample_name"))
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
  
  
  xlab("PC1 (75.9 %)") +
  ylab("PC2 (8.8 %)") +
  ggtitle("PCA_Proteome")



g1 <- g + scale_fill_manual(values = c("ALK" = "#147082", "ROS1" = "#8F0A2E"))

g1

library(FactoMineR)
library(factoextra)
gvar <- apply(mat, 1, var)
mostvargenes <- order(gvar, decreasing=TRUE)[1:nrow(mat)]
res_pca <- PCA(t(mat[mostvargenes,]), ncp=10, graph=TRUE)
fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main=("PCA_Proteome"))

nrow(mat)



#############################################################



###T-test---------------------------------------------

SamplePlan2 <- as.data.frame(t(Sample_Plan))
input_data <- rbind(SamplePlan2[3:4,],data_1pept2)

# Select samples
input_data1 <- as.data.frame(t(input_data))

input_data1A <-  subset(input_data1, input_data1[,"Oncofusion"] == "ALK" | input_data1[,"Oncofusion"] == "ROS1")

input_data2 <- Filter(function(x)!all(is.na(x)), input_data1A)

input_data2[,3:ncol(input_data2)] <- lapply(input_data2[,3:ncol(input_data2)], function(x) as.numeric(as.character(x)))
sapply(input_data2, class)

names(input_data2)[names(input_data2) == 'Oncofusion'] <- 'Proteome'




Sub <- "ALK"
TT <- input_data2[,2:ncol(input_data2)]

Tvalue <- sapply(colnames(TT)[2:ncol(TT)], function(x)
  tryCatch(
    t.test(TT[TT$Proteome==Sub, x], 
           TT[TT$Proteome!=Sub, x],
           paired = FALSE, var.equal = FALSE, alternative = "two.sided")$statistic,
    warning = function(w) return(NA),
    error = function (e) return(NA)
  ))

Tvalue <- as.data.frame(Tvalue)
rownames(Tvalue) <- colnames(TT)[2:ncol(TT)]

Pvalue <- sapply(colnames(TT)[2:ncol(TT)], function(x)
  tryCatch(
    t.test(TT[TT$Proteome==Sub, x], 
           TT[TT$Proteome!=Sub, x],
           paired = FALSE, var.equal = FALSE, alternative = "two.sided")$p.value,
    warning = function(w) return(NA),
    error = function (e) return(NA)
  ))

Pvalue <- as.data.frame(Pvalue)
rownames(Pvalue) <- colnames(TT)[2:ncol(TT)]


Meanvalue <- sapply(colnames(TT)[2:ncol(TT)], function(x)
  tryCatch(
    t.test(TT[TT$Proteome==Sub, x], 
           TT[TT$Proteome!=Sub, x],
           paired = FALSE, var.equal = FALSE, alternative = "two.sided")$estimate,
    warning = function(w) return(NA),
    error = function (e) return(NA)
  ))

Meanvalue <- as.data.frame(Meanvalue)
Meanvalue2 <- as.data.frame(t(Meanvalue))
colnames(Meanvalue2) <- c("Mean_Sub","Mean_Rest")
# Meanvalue2$Mean_ratio <- Meanvalue2$Mean_Sub/Meanvalue2$Mean_Rest
# Meanvalue2$Mean_Sub_unlog <- 10^(Meanvalue2$Mean_Sub) 
# Meanvalue2$Mean_Rest_unlog <- 10^(Meanvalue2$Mean_Rest) 
# Meanvalue2$FC <- Meanvalue2$Mean_Sub_unlog/Meanvalue2$Mean_Rest_unlog
Meanvalue2$LogFC <- Meanvalue2$Mean_Sub - Meanvalue2$Mean_Rest


Subgroup <- subset(TT, TT[,"Proteome"] == Sub)
Subgroup2 <- subset(Subgroup, select = -c(1))
Count_Subgroup <- as.data.frame(colSums(!is.na(Subgroup2)))
Tot_Subgroup <- nrow(Subgroup2)
Perc_Subgroup <- as.data.frame((Count_Subgroup*100)/Tot_Subgroup)
colnames(Perc_Subgroup) <- c("percentage_Sub")

Rest <- subset(TT, TT[,"Proteome"] != Sub)
Rest2 <- subset(Rest, select = -c(1))
Count_Rest <- as.data.frame(colSums(!is.na(Rest2)))
Tot_Rest <- nrow(Rest2)
Perc_Rest <- as.data.frame((Count_Rest*100)/Tot_Rest)
colnames(Perc_Rest) <- c("percentage_Rest")
Spec <- Perc_Subgroup/(Perc_Rest+Perc_Subgroup)
colnames(Spec) <- c("Specificity")

Protscore <- as.data.frame(cbind(Tvalue,Pvalue,Meanvalue2,Perc_Subgroup,Perc_Rest,Spec))
Rank_Unfilt_Protscore <- Protscore[order(Protscore$LogFC, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_ALK vs ROS1_without filters.txt",sep="\t",row.names=TRUE)



#Filter 
Protscore1 <- Protscore[(Protscore$percentage_Rest > 40),]
Protscore2 <- Protscore1[(Protscore1$percentage_Sub > 40),]
Protscore3 <-  Protscore2[(abs(Protscore2$LogFC) > 1),]
Rank_Protscore3 <- Protscore3[order(Protscore3$LogFC, decreasing = FALSE),, drop = FALSE]
Protscore4 <-  Protscore3[(Protscore3$Pvalue < 0.05),]
Rank_Protscore4 <- Protscore4[order(Protscore4$LogFC, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_Protscore4,"T-test_Oncofusion gliomas_ALK vs ROS1_LogFC_pval_filters.txt",sep="\t",row.names=TRUE)

Protscore8 <- Protscore[(Protscore$percentage_Sub > 60),]
Protscore9 <- Protscore8[(Protscore8$Specificity > 0.9),]
Rank_Protscore9 <- Protscore9[order(Protscore9$Specificity,Protscore9$percentage_Sub, decreasing = TRUE),, drop = FALSE]

Protscore10 <- Protscore[(Protscore$percentage_Rest > 60),]
Protscore11 <- Protscore10[(Protscore10$Specificity < 0.1),]
Rank_Protscore11 <- Protscore11[order(Protscore11$Specificity,Protscore11$percentage_Rest, decreasing = TRUE),, drop = FALSE]

Rank_Protscore12 <- rbind(Rank_Protscore9,Rank_Protscore11)
# write.table(Rank_Protscore12,"T-test_Oncofusion gliomas_ALK vs ROS1_Spec filter.txt",sep="\t",row.names=TRUE)



Rank_Protscore12["ROS1",]

#ENHANCED VOLCANO PLOTS----------------------------------------------------------------
library(EnhancedVolcano)


Prot <- row.names(Protscore2)
LogFC <- Protscore2$LogFC
Pvalue <- Protscore2$Pvalue
logPvalue <- -log10(Pvalue)

df <- as.data.frame(cbind(Prot,LogFC,Pvalue,logPvalue))
row.names(df) <- df[,"Prot"]
sapply(df, class)

df[,2:ncol(df)] <- lapply(df[,2:ncol(df)], function(x) as.numeric(as.character(x)))
sapply(df, class)

df[is.na(df)] <- 1

IA2T2 <- df[(df$LogFC > 2.5),]
IA2T1 <- df[(df$LogFC < -2.1),]
IA2T1 <- IA2T1[(IA2T1$Pvalue < 0.001),]
IA2T2 <- IA2T2[(IA2T2$Pvalue < 0.001),]
IA2T0 <- df[(df$Pvalue < 10e-50),]
IA2T0 <- IA2T0[(abs(IA2T0$LogFC) > 1),]


EnhancedVolcano(df,
                lab = rownames(df),
                x = 'LogFC',
                y = 'Pvalue',
                ylim = c(0,60),
                xlim = c(-4,4),
                selectLab = c(rownames(IA2T2),rownames(IA2T1),rownames(IA2T0),"ALK","ROS1"),
                title = 'ALK vs ROS1 (Oncofusion gliomas)',
                pCutoff = 0.001,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4,
                col=c('grey', 'grey', 'grey', 'red4'),
                colAlpha = 1,
                drawConnectors = TRUE,
                maxoverlapsConnectors = Inf,
                widthConnectors = 0.75,
                colConnectors = 'black',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE
)


dev.off()




###Boxplot----------------------

set.seed(123)
Protein <- "ROS1"
data_Protein <- as.data.frame(input_data2[,c("Proteome",Protein)])
colnames(data_Protein) <- c("Proteome", "Protein")

data_ProteinB <- data_Protein
data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
                                levels = c("ALK", "ROS1"),
                                labels = c("ALK", "ROS1"))


par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
        lwd = 2, xlab = 'Oncofusion', ylab = 'Protein level', col = c("blue","red"),
        ylim=c(-3.2,3.2),
        outline=FALSE)

stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("darkblue","darkred"))



dev.off()

#Protein selection by highest SD-------------------------------------------------------------------------------------------
library(multiClust)
#Probe Ranking
mat_500 <- probe_ranking(input=NULL, probe_number=500,
                         probe_num_selection="Fixed_Probe_Num", data.exp=mat, method="SD_Rank")
library(genefilter)
#Select proteins with highest IQR
mat_300B <- varFilter(mat, var.func=IQR, var.cutoff=0.887, filterByQuantile=TRUE)

#Select proteins with highest SD
vars <- apply(mat, 1, sd)
mat500B <- mat[vars > quantile(vars, 0.8), ] 

# sd(mat300C["AKAP8",])
# IQR(mat300C["AKAP8",])

mat <- as.matrix(mat500B)


### Most variable proteins by condition / group
# mat <- prot.10na
mat <- mean.centered.prot.10na

Pval <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Sample_name)
  return(summary(trans)[[1]][1,5])
}

t.mat <- apply(mat, 1, Pval)
t.mat <- t.mat[t.mat <= 1e-26]
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
Nclust <- 2
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
              column_title = "HC TOP1110 most variable proteins",
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
    proteomic <- cbind(clu, paste("cluster_PROT_", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(row.names(mat[hc[[i]],])))
    clu <- cbind(clu, paste("cluster_PROT_", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}


proteomic
colnames(proteomic) <- c("Gene","Cluster_PROT") 
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper")
# write.table(proteomic,"list_Oncofusion_cluster_PROT_TOP1110.txt",sep="\t",row.names=FALSE)




# ComplexHeatmap for specific proteins-----------------------------


# PTM_SPEC <- read.delim("T-test_dTAG vs T0_Spec filter.txt")
PTM_SPEC <- Rank_Protscore12
PTM_SPEC3 <-  PTM_SPEC


list_SPEC <- row.names(PTM_SPEC3)

data_1pept2_sub <- t(input_data2)

data_1pept_SPEC <- as.data.frame(data_1pept2_sub[row.names(data_1pept2_sub) %in% list_SPEC , ])
data_1pept_SPEC[,] <- sapply(data_1pept_SPEC[,], function(x) as.numeric(as.character(x)))
sapply(data_1pept_SPEC, class)

boxplot(data_1pept_SPEC)
data_1pept_SPEC1 <- data_1pept_SPEC[match(list_SPEC, row.names(data_1pept_SPEC)),]
mat <- as.matrix(data_1pept_SPEC1)
# write.table(mat,"Specific_proteins.txt",sep="\t",row.names=TRUE)

SamplePlan2 <- as.data.frame(t(Sample_Plan))
HC <- as.data.frame(mat)
SamplePlan2B <- SamplePlan2[names(HC)]
Condition <- SamplePlan2B["Oncofusion",]
mat_G1 <- as.matrix(Condition)
Cell_Line <- SamplePlan2B["Sample_name",]
mat_CL <- as.matrix(Cell_Line)


ht <- Heatmap(mat, na_col = "lightgrey", 
              # col = colorRamp2(c(-2, 1, 0, 1, 2), c("darkblue","blue", "white", "red","darkred")),
              column_dend_reorder = FALSE,
              column_title = "Specific proteins_ALK vs ROS1",
              name = "Protein level", 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              # column_split = Nclust,
              show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 10), 
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


dev.off()

# Monte Carlo Clustering-------------------------------
library(M3C)
library(NMF) # loading for aheatmap plotting function
library(gplots) # loading this for nice colour scale
library(ggsci) # more cool colours

pca(mat,legendtextsize = 10,axistextsize = 10,dotsize=2)

res <- M3C(mat, maxK = 10, method = 1, removeplots = FALSE, iters=100, cores=1, seed = 456, clusteralg = "hc")
# res$scores
# res$plots
# M3C

Nclust <- 7
M3C <- res[["realdataresults"]][[Nclust]][["consensus_matrix"]]
data <- res$realdataresults[[Nclust]]$ordered_data # this is the data
annon <- res$realdataresults[[Nclust]]$ordered_annotation
# pca(data,labels=annon$consensuscluster,legendtextsize = 10,axistextsize = 10,dotsize=2)
colnames(M3C) <- colnames(data)
mat_M3C <- as.matrix(M3C)

SamplePlan2 <- as.data.frame(t(Sample_Plan))
M3C <- as.data.frame(M3C)
SamplePlan2 <- SamplePlan2[names(M3C)]
Group <- SamplePlan2["Oncofusion",]
mat_G1 <- as.matrix(Group)



Prot_M3C <- res[["realdataresults"]][[Nclust]][["ordered_annotation"]]
str(Prot_M3C)
Prot_M3C[] <- lapply(Prot_M3C, function(x) as.character(as.factor(x)))
sapply(Prot_M3C,class)
Prot_M3C$New_name <- row.names(Prot_M3C)


# count_prot <- prot.10na
# count_prot1 <- count_prot
# count_prot2 <- as.data.frame(colSums(!is.na(count_prot1)))
# names(count_prot2)[names(count_prot2) == 'colSums(!is.na(count_prot1))'] <- 'Total_Proteins'
# str(count_prot2)
# count_prot3 <- as.data.frame(t(count_prot2))
# M3C <- as.data.frame(M3C)
# count_prot3 <- count_prot3[names(M3C)]
# mat_count_Prot <- as.matrix(count_prot3)


ht <- Heatmap(mat_M3C, na_col = "white", 
              col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue", "lightblue", "beige", "red1", "darkred")),
              column_dend_reorder = FALSE,
              column_title = "Unsupervised M3C clustering_6188 proteins (k=2)",
              name = "Consensus score", 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              # column_split = 8,
              # show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 5), 
              column_names_gp = gpar(fontsize = 12), 
              top_annotation = HeatmapAnnotation(
                Oncofusion = mat_G1[1,],
                col = list(
                  Oncofusion = c("ALK" = "blue", "ROS1" = "red")
                )
              ))

ht = draw(ht)


proteomic5C <- SamplePlan2[c(3,4),]




# PCA with subgroups------------------------------------

proteomic9 <- as.data.frame(rbind(proteomic5C,mat[, colnames(proteomic5C)]))

proteomic10 <- as.data.frame(t(proteomic9))

colnames(proteomic10[,1:4])

mat_PCA <- proteomic10 %>% select(-contains("Sample_name")) %>% select(-contains("Oncofusion"))
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
# mat_PCA3 <- cbind(mat_PCA2,proteomic10$Oncofusion)
mat_PCA3 <- cbind(mat_PCA2,proteomic10$Sample_name)



my_data_PCA <- cbind(mat_PCA3, myPr$x)

my_data_PCA2 <- as.data.frame(my_data_PCA)

# names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Oncofusion'] <- 'Oncofusion'
names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Sample_name'] <- 'Sample'
# names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Condition'] <- 'Condition'


g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2))
                              , col = Cell_Line, fill = Cell_Line
                              )) +
  # stat_ellipse(geom = "polygon", col = "black", alpha = 0.3) +
  geom_point(size = 5, shape = 21, col = "black") +
  theme_light(base_size=20) +
  
  
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA_6188 proteins")
  
g

g1 <- g + scale_fill_manual(values = c("ALK" = "blue", "ROS1" = "red"))

g1

# + geom_text(col = "black", size = 4, check_overlap = FALSE, vjust = 0, nudge_y = 0.3)



#PCA in 3D (FORMAT1)--------------------------------------------------------------------------------------
library(rgl)
library(car)


PC1 <- my_data_PCA2[,"PC1"]
PC2 <- my_data_PCA2[,"PC2"]
PC3 <- my_data_PCA2[,"PC3"]
Cluster <- as.factor(my_data_PCA2[,"Cell_Line"])

scatter3d(x = PC1, y = PC2, z = PC3, groups = Cluster,
          surface=FALSE, grid = FALSE, ellipsoid = FALSE,
          # surface.col = c("ALK" = "blue", "ROS1" = "red"),
          axis.col = c("black", "black", "black", "black", "black", "black", "black", "black", "black"),sphere.size=1.15, radius=1)


rgl.snapshot(filename = "PCA_3D_oncofusion.png")

#PCA in 3D (FORMAT2)--------------------------------------------------------------------------------------
library(plot3D)

colors <- c("green", "blue")
colors <- colors[as.numeric(my_data_PCA2$Proteomic)]
p3d <- plot3d(my_data_PCA2$PC1, my_data_PCA2$PC2, my_data_PCA2$PC3, col= colors, size=1, type='s' 
)
text3d(my_data_PCA2$PC1, my_data_PCA2$PC2, my_data_PCA2$PC3+1,
       texts=c(rownames(my_data_PCA2)),col= colors , cex= 0.7, pos=3)

rgl.snapshot(filename = "PCA3Dbis.png")





#UMAP-----------------------------------------------------------------------------------

library(umap)
set.seed(456)
umap <- umap(mat_PCA2,n_neighbors = 5L)


df <- data.frame(x = umap$layout[,1],
                 y = umap$layout[,2],
                 Group = PDX)





t <- ggplot(df, aes(x, y, colour = Group)) +
  geom_point(size = 6,pch = 16) +
  theme_light(base_size=20) +
  
  
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("UMAP_Proteome_Osteosarcoma") 


# t1 <- t + scale_colour_manual(values = c("T0" = "green", "T4h.DMSO" = "pink", "T4h.dTAG" = "red","T8h.DMSO" = "cyan","T8h.dTAG" = "blue"))
  
t



delete <- umap[["data"]]
