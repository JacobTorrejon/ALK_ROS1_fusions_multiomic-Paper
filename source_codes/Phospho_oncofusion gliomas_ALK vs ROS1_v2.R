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


#OPEN DATA FILE----------------------------------------------------------------------------------------------
# Open matrix
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024/phosphoMS ROS1 and ALK SerThr and Tyr/ALK")
matrix1 <- read.delim("onlyLocalizedSites_o33970_TMTPhospho_STplusY_towardsMOFA_PepQuantMat.txt", header = TRUE)
rownames(matrix1) = make.names(matrix1$GeneNameNphosphoSite, unique = TRUE)

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024/phosphoMS ROS1 and ALK SerThr and Tyr/ROS1")
matrix2 <- read.delim("onlyLocalizedSites_o31555_TMTPhospho_STplusY_towardsMOFA_PepQuantMat.txt", header = TRUE)
rownames(matrix2) = make.names(matrix2$GeneNameNphosphoSite, unique = TRUE)


matrix <- merge(x= matrix1[,-c(1:4,21:25)], y=matrix2[,-c(1:4,23:27)], by='row.names', all = TRUE)

list_prot <- matrix$Row.names
row.names(matrix) <- matrix$Row.names
boxplot(matrix[,2:35])

data_1pept <- matrix[,-c(1)]
data_1pept1 <- as.data.frame(scale(data_1pept, center = TRUE, scale = TRUE))
boxplot(data_1pept1)


#Checking %NAs
sum(is.na(data_1pept1))
sum(!is.na(data_1pept1))
sum(is.na(data_1pept1))*100/(ncol(data_1pept1)*nrow(data_1pept1))

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024")
Sample_Plan <- read.csv("Sample_Plan_ALKvsROS1.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name

Sample_ID <- Sample_Plan$Sample_ID

data_1pept2 <- data_1pept1[,Sample_ID]
colnames(data_1pept2) <- row.names(Sample_Plan)
sapply(data_1pept2, class)

par(mar = c(12, 4, 4, 2) + 0.1);
boxplot(data_1pept2,las=2, ylab = "PTM level")

dev.off()

#DIAGNOSIS PLOTS FOR PCA IMPUTATION---------------------------------------------------------------------------------------------------

#Load data
prot <- as.matrix(data_1pept2)  ## Proteins in rows, samples in columns 
nzero.prot <- apply(prot, MARGIN = 1, FUN = function(x) length(x[is.na(x)]) )

# Plot N.missing values vs N.quantified proteins
df.nzero.prot <- data.frame(table(nzero.prot))
x.df.nzero.prot <- as.numeric(as.character(df.nzero.prot[,1]))
x.df.nzero.prot2 <- x.df.nzero.prot*100/ncol(prot)

y.df.nzero.prot <- cumsum(df.nzero.prot[,2])
plot(x.df.nzero.prot2, y.df.nzero.prot, pch = 16,  main = "", xlab = "% Missing values", 
     ylab = "Quantified proteins", cex=2, ylim = c(0,11000) )
grid()

## Create a matrix with at most 34 NAs per protein
na.th = ncol(prot)*0.01
prot.10na <- prot[which(nzero.prot<=na.th),] 
boxplot(prot.10na)
mean.centered.prot.10na <- scale(prot.10na, center = TRUE, scale = TRUE)

boxplot(mean.centered.prot.10na)


# saveRDS(prot.10na,"Phospho matrix_ROS1_Normalized.rds")

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
#   ylim(0,16000)+
#   facet_wrap(~ var, scales = "free") +
#   theme_bw()




# mat <- t(prot.10na_imput)
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
  
  
  xlab("PC1 (74.9 %)") +
  ylab("PC2 (9.2 %)") +
  ggtitle("PCA_Phospho")



g1 <- g + scale_fill_manual(values = c("ALK" = "#147082", "ROS1" = "#8F0A2E"))

g1

library(FactoMineR)
library(factoextra)
gvar <- apply(mat, 1, var)
mostvargenes <- order(gvar, decreasing=TRUE)[1:nrow(mat)]
res_pca <- PCA(t(mat[mostvargenes,]), ncp=10, graph=TRUE)
fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main=("PCA_Phospho"))

nrow(mat)

##############################################################################################################
# 8. KINASE ACTIVITY FROM MATRIX                                                                                   #
##############################################################################################################

library(OmnipathR)
library(decoupleR)

# Pval <- function(x)
# {
#   trans <- aov(x ~ Sample_Plan$Sample_name)
#   return(summary(trans)[[1]][1,5])
# }
# 
# 
# phospho_p <- prot.10na
# t.phospho <- apply(phospho_p, 1, Pval)
# t.phospho <- t.phospho[t.phospho <= 0.01]
# phospho_p <- phospho_p[names(t.phospho),]


Pval2 <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Kinase_Dead)
  return(summary(trans)[[1]][1,5])
}


phospho_p <- prot.10na
t.phospho <- apply(phospho_p, 1, Pval2)
t.phospho <- t.phospho[t.phospho <= 0.01]
phospho_p <- phospho_p[names(t.phospho),]




phospho_p <- as.data.frame(phospho_p)
phospho_p2 <-tibble::rownames_to_column(phospho_p, "PTM")



MOFA_F1B <- phospho_p2
# MOFA_F1B <- Protscore4[,c("PTM","LogFC")]
# colnames(MOFA_F1B) <- c("PTM","Tvalue")



MOFA_F1C <- MOFA_F1B %>% separate(PTM, c("Gene","site1", "site2", "site3"),fill = "right")
# MOFA_F1C["H1.6.S181","Gene"] <- c("H1-6")
# MOFA_F1C["H1.6.S181","site1"] <- c("S181")
# MOFA_F1C["H1.6.S181","site2"] <- NA


MOFA_F1C$PTM <- with(MOFA_F1C, paste(Gene,site1, sep="_"))
list_Samples <- c("GOPC.ROS1_1","GOPC.ROS1_2","GOPC.ROS1_3","GOPC.ROS1_KD_1","GOPC.ROS1_KD_2","CLIP1.ROS1_1","GOPC.ROS1_KD_3",  
                  "CLIP1.ROS1_2","CLIP1.ROS1_3","CLIP1.ROS1_KD_1","CLIP1.ROS1_KD_2","KIF21A.ROS1_1","CLIP1.ROS1_KD_3","KIF21A.ROS1_2",
                  "KIF21A.ROS1_3","KIF21A.ROS1_KD_1","KIF21A.ROS1_KD_2","KIF21A.ROS1_KD_3")

MOFA_PTM1 <- MOFA_F1C[,c("PTM",list_Samples)]

MOFA_PTM2 <- MOFA_F1C[,c("Gene","site2",list_Samples)] %>% drop_na(site2)
MOFA_PTM2$PTM <- with(MOFA_PTM2, paste(Gene,site2, sep="_"))
MOFA_PTM2 <- MOFA_PTM2[,c("PTM",list_Samples)]

MOFA_PTM3 <- MOFA_F1C[,c("Gene","site3",list_Samples)] %>% drop_na(site3)
MOFA_PTM3$PTM <- with(MOFA_PTM3, paste(Gene,site3, sep="_"))
MOFA_PTM3 <- MOFA_PTM3[,c("PTM",list_Samples)]

# MOFA_PTM4 <- MOFA_F1C[,c("Gene","site4","Tvalue")] %>% drop_na(site4)
# MOFA_PTM4$PTM <- with(MOFA_PTM4, paste(Gene,site4, sep="_"))
# MOFA_PTM4 <- MOFA_PTM4[,c("PTM","Tvalue")]

MOFA_PTM <- rbind(MOFA_PTM1,MOFA_PTM2,MOFA_PTM3)
MOFA_PTM$Mean <- rowMeans(MOFA_PTM[,c(-1)])

MOFA_PTM_Ranked <- MOFA_PTM[order(MOFA_PTM$Mean, decreasing = TRUE), ]
MOFA_PTM_Ranked_Nodup <- MOFA_PTM_Ranked[!duplicated(MOFA_PTM_Ranked$PTM), ]
row.names(MOFA_PTM_Ranked_Nodup) <- MOFA_PTM_Ranked_Nodup$PTM


counts <- as.matrix(MOFA_PTM_Ranked_Nodup[,list_Samples])
head(counts)


design <- Sample_Plan[,c("Replica_name","Sample_name","Kinase_Dead","Oncofusion")]
colnames(design) <- c("sample","Cell_Line","Kinase_Dead","Oncofusion")
design


### Next, we can load the prior knowledge interactions, composed by kinase-target relationships
uniprot_kinases <- OmnipathR::import_omnipath_annotations(resources = "UniProt_keyword") %>%
  dplyr::filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
  distinct() %>%
  pull(genesymbol) %>%
  unique()
omnipath_ptm <- OmnipathR::get_signed_ptms() %>%
  dplyr::filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
  dplyr::filter(!(stringr::str_detect(sources, "ProtMapper") & n_resources == 1)) %>%
  dplyr::mutate(p_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset),
                mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
  dplyr::transmute(p_site, enzyme_genesymbol, mor) %>%
  dplyr::filter(enzyme_genesymbol %in% uniprot_kinases)

omnipath_ptm$likelihood <- 1


#we remove ambiguous modes of regulations
omnipath_ptm$id <- paste(omnipath_ptm$p_site,omnipath_ptm$enzyme_genesymbol, sep ="")
omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id),]
omnipath_ptm <- omnipath_ptm[,-5]

#On a final step, we run viper to get the Kinase activities from the phosphoproteomic data
names(omnipath_ptm)[c(1,2)] <- c("target","tf")

sample_acts <- run_ulm(mat=counts, net=omnipath_ptm, .source='tf', .target='target',
                       .mor='mor', minsize = 5)


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
t.sample_acts <- t.sample_acts[t.sample_acts <= 5e-2]
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


Nclust <- 2
# Plot
pheatmap(sample_acts_mat2, border_color = NA, color=my_color, breaks = my_breaks,
         cluster_cols=F, cluster_rows=T, 
         scale="row",
         clustering_distance_cols="euclidean", clustering_method="ward.D",
         # cutree_cols = Nclust,
         cutree_rows = Nclust,
         annotation_col=data.frame(metadata),
         annotation_colors = annoCol,
         show_rownames= TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         fontsize = 8,
         main = "Kinase activity_Most variables KD vs CTRL(ALL)") 


## Boxplot
input_data2 <- cbind(design,sample_acts_mat)
set.seed(123)
TF <- "MTOR"
data_Protein <- as.data.frame(input_data2[,c("Cell_Line",TF)])
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


#######################################################
###T-test---------------------------------------------
#######################################################

SamplePlan2 <- as.data.frame(t(Sample_Plan))
input_data <- rbind(SamplePlan2[c(4,5),],data_1pept2)
input_data1 <- as.data.frame(t(input_data))

input_data2 <- Filter(function(x)!all(is.na(x)), input_data1)

input_data2[,3:ncol(input_data2)] <- lapply(input_data2[,3:ncol(input_data2)], function(x) as.numeric(as.character(x)))
sapply(input_data2, class)

names(input_data2)[names(input_data2) == 'Oncofusion'] <- 'Proteome'
# names(input_data2)[names(input_data2) == 'KD'] <- 'Proteome'





####---------------------------------------------------------------------


Sub <- "ALK"


TT <- input_data2[,c(2:ncol(input_data2))]

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
Protscore$PTM <- row.names(Protscore)
Protscore <- Protscore[,c(9,1:8)]
Rank_Unfilt_Protscore <- Protscore[order(Protscore$LogFC, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_GOPC.ROS1 vs CLIP1.ROS1_without filters.txt",sep="\t",row.names=FALSE)
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_GOPC.ROS1 vs KIF21A.ROS1_without filters.txt",sep="\t",row.names=FALSE)
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_CLIP1.ROS1 vs KIF21A.ROS1_without filters.txt",sep="\t",row.names=FALSE)
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_GOPC.ROS1 vs GOPC.ROS1_KD_without filters.txt",sep="\t",row.names=FALSE)
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_CLIP1.ROS1 vs CLIP1.ROS1_KD_without filters.txt",sep="\t",row.names=FALSE)
# write.table(Rank_Unfilt_Protscore,"T-test_Oncofusion gliomas_KIF21A.ROS1 vs KIF21A.ROS1_KD_without filters.txt",sep="\t",row.names=FALSE)



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
# write.table(Rank_Protscore12,"Specific_PTMs_Oncofusion gliomas_ALK vs ROS1_Spec filter.txt",sep="\t",row.names=TRUE)



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

IA2T2 <- df[(df$LogFC > 5),]
IA2T1 <- df[(df$LogFC < -3),]
IA2T1 <- IA2T1[(IA2T1$Pvalue < 0.05),]
IA2T2 <- IA2T2[(IA2T2$Pvalue < 0.05),]
IA2T0 <- df[(df$Pvalue < 2e-40),]
IA2T0 <- IA2T0[(abs(IA2T0$LogFC) > 0.9),]


EnhancedVolcano(df,
                lab = rownames(df),
                x = 'LogFC',
                y = 'Pvalue',
                ylim = c(0,50),
                xlim = c(-6,6),
                selectLab = c(rownames(IA2T2),rownames(IA2T1),rownames(IA2T0),"ALK","ROS1"),
                title = 'KIF21A.ROS1 vs KIF21A.ROS1_KD (Oncofusion gliomas)',
                pCutoff = 0.05,
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
# write.table(mat,"Specific_proteins.txt",sep="\t",row.names=TRUE)

SamplePlan_SPEC <- Sample_Plan %>% arrange(Oncofusion,Sample_name)
t <- row.names(SamplePlan_SPEC)
SamplePlan2B <- as.data.frame(t(SamplePlan_SPEC))
Condition <- SamplePlan2B["Oncofusion",]
mat_G1 <- as.matrix(Condition)
Cell_Line <- SamplePlan2B["Sample_name",]
mat_CL <- as.matrix(Cell_Line)
mat <- as.matrix(data_1pept_SPEC1[,t])


ht <- Heatmap(mat, na_col = "lightgrey", 
              # col = colorRamp2(c(-2, 1, 0, 1, 2), c("darkblue","blue", "white", "red","darkred")),
              column_dend_reorder = FALSE,
              column_title = "Specific PTMs_ALK vs ROS1",
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





###Boxplot----------------------

set.seed(123)
# Protein <- "ROS1.Y2110.Y2114"
Protein <- "ROS1.Y2274"
data_Protein <- as.data.frame(input_data2[,c("Proteome",Protein)])
colnames(data_Protein) <- c("Proteome", "Protein")

data_ProteinB <- data_Protein
# data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
#                                  levels = c("GOPC.ROS1", "GOPC.ROS1_KD","CLIP1.ROS1", "CLIP1.ROS1_KD","KIF21A.ROS1", "KIF21A.ROS1_KD"),
#                                  labels = c("GOPC.ROS1", "GOPC.ROS1_KD","CLIP1.ROS1", "CLIP1.ROS1_KD","KIF21A.ROS1", "KIF21A.ROS1_KD"))

par(mar = c(12, 4, 4, 2) + 0.6);
par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
        lwd = 2, xlab = '', ylab = 'PTM level', col = c("red","pink","yellow2","yellow1","orange2","orange1"),las=2,
        # ylim=c(6,16),
        outline=FALSE)

stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=8, pch = 20, col = c("black","black","black","black","black","black"))


library(ggpubr)
df.factor2 <- data_ProteinB
colnames(df.factor2) <- c("Cell_Line","value")

p.bx.subtprot.G3B <- ggplot(df.factor2, aes(x=Cell_Line, y=value)) +  
  # geom_violin(aes(fill=Cell_Line), alpha = 1, lwd = 1.2, show.legend = FALSE) +
  geom_boxplot(aes(group=Cell_Line, fill = Cell_Line),
               width=0.4, show.legend = FALSE, alpha = 0.5,
               # position=position_dodge(width=0.35), pch=21, size=1,  stroke = 1.5
  ) + 
  # geom_point(colour = "#4D4D4D", size = 3)+
  geom_point(aes( fill = Cell_Line), size = 3,pch = 21, show.legend = FALSE) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) + 
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20), axis.text = element_text(size = 20, angle = 90)) + ylim(10,16.8)+
  geom_signif(comparisons = list(c("CLIP1:ROS1", "GOPC:ROS1"),c("CLIP1:ROS1", "KIF21A:ROS1")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = 0.014,
              step_increase = 0.05,
              extend_line = 0,
              tip_length = 0.01,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "PTM level")) + ggtitle("ROS1.Y2274")





# data_ProteinB <- data_Protein
# data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
#                                  levels = c("KIF21A.ROS1", "KIF21A.ROS1_KD"),
#                                  labels = c("KIF21A.ROS1", "KIF21A.ROS1_KD"))
# 
# par(mar = c(12, 4, 4, 2) + 0.6);
# par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
# boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
#         lwd = 2, xlab = '', ylab = 'PTM level', col = c("orange2","orange1"),las=2,
#         # ylim=c(-3.2,3.2),
#         outline=FALSE)
# 
# stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
#            method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("orange4","orange2"))



data_ProteinB <- data_Protein
data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
                                 levels = c("CLIP1.ROS1", "CLIP1.ROS1_KD"),
                                 labels = c("CLIP1.ROS1", "CLIP1.ROS1_KD"))

par(mar = c(12, 4, 4, 2) + 0.6);
par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
        lwd = 2, xlab = '', ylab = 'PTM level', col = c("yellow2","yellow1"),las=2,
        # ylim=c(-3.2,3.2),
        outline=FALSE)

stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("yellow4","yellow2"))



data_ProteinB <- data_Protein
data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
                                levels = c("GOPC.ROS1", "GOPC.ROS1_KD"),
                                labels = c("GOPC.ROS1", "GOPC.ROS1_KD"))

par(mar = c(12, 4, 4, 2) + 0.6);
par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
        lwd = 2, xlab = '', ylab = 'PTM level', col = c("red","pink"),las=2,
        # ylim=c(-3.2,3.2),
        outline=FALSE)

stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("darkred","tomato"))


# data_ProteinB <- data_Protein
# data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
#                                 levels = c("GOPC.ROS1", "CLIP1.ROS1"),
#                                 labels = c("GOPC.ROS1", "CLIP1.ROS1"))
# 
# par(mar = c(12, 4, 4, 2) + 0.6);
# par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
# boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
#         lwd = 2, xlab = '', ylab = 'PTM level', col = c("red","yellow2"),las=2,
#         # ylim=c(-3.2,3.2),
#         outline=FALSE)
# 
# stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
#            method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("darkred","yellow4"))


# data_ProteinB <- data_Protein
# data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
#                                 levels = c("GOPC.ROS1", "KIF21A.ROS1"),
#                                 labels = c("GOPC.ROS1", "KIF21A.ROS1"))
# 
# par(mar = c(12, 4, 4, 2) + 0.6);
# par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
# boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
#         lwd = 2, xlab = '', ylab = 'PTM level', col = c("red","orange2"),las=2,
#         # ylim=c(-3.2,3.2),
#         outline=FALSE)
# 
# stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
#            method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("darkred","orange3"))



data_ProteinB <- data_Protein
data_ProteinB$Proteome <- factor(data_ProteinB$Proteome, exclude = NULL,
                                 levels = c("CLIP1.ROS1", "KIF21A.ROS1"),
                                 labels = c("CLIP1.ROS1", "KIF21A.ROS1"))

par(mar = c(12, 4, 4, 2) + 0.6);
par(cex.lab=1.5, cex.axis=1.4, cex.main = 2.0)
boxplot(Protein ~ Proteome, data = data_ProteinB, main = Protein, boxwex = 0.5,
        lwd = 2, xlab = '', ylab = 'PTM level', col = c("yellow2","orange2"),las=2,
        # ylim=c(-3.2,3.2),
        outline=FALSE)

stripchart(Protein ~ Proteome, vertical = TRUE, data = data_ProteinB,
           method = "jitter", add = TRUE, lwd=15, pch = 20, col = c("yellow4","orange3"))


dev.off()


##############################################################################################################
# 8. KINASE ACTIVITY                                                                                         #
##############################################################################################################

library(OmnipathR)
library(decoupleR)

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024/phosphoMS ROS1 and ALK SerThr and Tyr/ROS1")
# Protscore <- read.delim("T-test_Oncofusion gliomas_GOPC.ROS1 vs CLIP1.ROS1_without filters.txt")
# Protscore <- read.delim("T-test_Oncofusion gliomas_GOPC.ROS1 vs KIF21A.ROS1_without filters.txt")
# Protscore <- read.delim("T-test_Oncofusion gliomas_CLIP1.ROS1 vs KIF21A.ROS1_without filters.txt")
# Protscore <- read.delim("T-test_Oncofusion gliomas_GOPC.ROS1 vs GOPC.ROS1_KD_without filters.txt")
Protscore <- read.delim("T-test_Oncofusion gliomas_CLIP1.ROS1 vs CLIP1.ROS1_KD_without filters.txt")
# Protscore <- read.delim("T-test_Oncofusion gliomas_KIF21A.ROS1 vs KIF21A.ROS1_KD_without filters.txt")
row.names(Protscore) <- Protscore$PTM

Protscore1 <- Protscore[(Protscore$percentage_Rest > 40),]
Protscore2 <- Protscore1[(Protscore1$percentage_Sub > 40),]
Protscore3 <-  Protscore2[(abs(Protscore2$LogFC) > 0.5),]
Rank_Protscore3 <- Protscore3[order(Protscore3$LogFC, decreasing = FALSE),, drop = FALSE]
Protscore4 <-  Protscore3[(Protscore3$Pvalue < 0.05),]

MOFA_F1B <- Protscore4[,c("PTM","LogFC")]
colnames(MOFA_F1B) <- c("PTM","Tvalue")

# setwd("~/Documents/Curie_Jacob/Figures_MB multiomic paper/MOFA_RNA_Prot_Phospho_Fig3/G3")
# MOFA_F1 <- read.delim("Rank_weights_Phospho_Factor1.txt")
# phospho_input <- MOFA_F1 %>%
#   rownames_to_column(var="Protein")
# colnames(phospho_input) <- c("Protein","LogFC")
# row.names(phospho_input) <- phospho_input$Protein
# phospho_input2 <- phospho_input[(abs(phospho_input$LogFC) > 0.2),]



MOFA_F1C <- MOFA_F1B %>% separate(PTM, c("Gene","site1", "site2", "site3"),fill = "right")
# MOFA_F1C["H1.6.S181","Gene"] <- c("H1-6")
# MOFA_F1C["H1.6.S181","site1"] <- c("S181")
# MOFA_F1C["H1.6.S181","site2"] <- NA


MOFA_F1C$PTM <- with(MOFA_F1C, paste(Gene,site1, sep="_"))
MOFA_PTM1 <- MOFA_F1C[,c("PTM","Tvalue")]

MOFA_PTM2 <- MOFA_F1C[,c("Gene","site2","Tvalue")] %>% drop_na(site2)
MOFA_PTM2$PTM <- with(MOFA_PTM2, paste(Gene,site2, sep="_"))
MOFA_PTM2 <- MOFA_PTM2[,c("PTM","Tvalue")]

MOFA_PTM3 <- MOFA_F1C[,c("Gene","site3","Tvalue")] %>% drop_na(site3)
MOFA_PTM3$PTM <- with(MOFA_PTM3, paste(Gene,site3, sep="_"))
MOFA_PTM3 <- MOFA_PTM3[,c("PTM","Tvalue")]

# MOFA_PTM4 <- MOFA_F1C[,c("Gene","site4","Tvalue")] %>% drop_na(site4)
# MOFA_PTM4$PTM <- with(MOFA_PTM4, paste(Gene,site4, sep="_"))
# MOFA_PTM4 <- MOFA_PTM4[,c("PTM","Tvalue")]

MOFA_PTM <- rbind(MOFA_PTM1,MOFA_PTM2,MOFA_PTM3)
MOFA_PTM_Ranked <- MOFA_PTM[order(MOFA_PTM$Tvalue, decreasing = TRUE), ]
MOFA_PTM_Ranked_Nodup <- MOFA_PTM_Ranked[!duplicated(MOFA_PTM_Ranked$PTM), ]
row.names(MOFA_PTM_Ranked_Nodup) <- MOFA_PTM_Ranked_Nodup$PTM


setwd("~/Documents/Curie_Jacob/Scripts to explore/kinase_tf_mini_tuto-main")
library(here)
source(here("code/utils.R"))

phospho_differential_analysis <- as.data.frame(MOFA_PTM_Ranked_Nodup[,c(2)])
row.names(phospho_differential_analysis) <- MOFA_PTM_Ranked_Nodup$PTM
colnames(phospho_differential_analysis) <- c("MOFA_score")

plot_top_features(phospho_differential_analysis, n_top = 40) +
  ggtitle('Phosphosite space')


### Next, we can load the prior knowledge interactions, composed by kinase-target relationships
uniprot_kinases <- OmnipathR::import_omnipath_annotations(resources = "UniProt_keyword") %>%
  dplyr::filter(value == "Kinase" & !grepl("COMPLEX", uniprot)) %>%
  distinct() %>%
  pull(genesymbol) %>%
  unique()
omnipath_ptm <- OmnipathR::get_signed_ptms() %>%
  dplyr::filter(modification %in% c("dephosphorylation","phosphorylation")) %>%
  dplyr::filter(!(stringr::str_detect(sources, "ProtMapper") & n_resources == 1)) %>%
  dplyr::mutate(p_site = paste0(substrate_genesymbol, "_", residue_type, residue_offset),
                mor = ifelse(modification == "phosphorylation", 1, -1)) %>%
  dplyr::transmute(p_site, enzyme_genesymbol, mor) %>%
  dplyr::filter(enzyme_genesymbol %in% uniprot_kinases)

omnipath_ptm$likelihood <- 1


#we remove ambiguous modes of regulations
omnipath_ptm$id <- paste(omnipath_ptm$p_site,omnipath_ptm$enzyme_genesymbol, sep ="")
omnipath_ptm <- omnipath_ptm[!duplicated(omnipath_ptm$id),]
omnipath_ptm <- omnipath_ptm[,-5]

#On a final step, we run viper to get the Kinase activities from the phosphoproteomic data
names(omnipath_ptm)[c(1,2)] <- c("target","tf")

kin_activity <- run_ulm(mat=phospho_differential_analysis, net=omnipath_ptm, .source='tf', .target='target',
                        .mor='mor', minsize = 5)


##### PLOT TOP TF
TF <- kin_activity

n_tfs <- 21

f_TF <- TF %>%
  filter(statistic == 'ulm') %>%
  mutate(rnk = NA)

# Filter top TFs in both signs
msk <- f_TF$score > 0
f_TF[msk, 'rnk'] <- rank(-f_TF[msk, 'score'])
f_TF[!msk, 'rnk'] <- rank(-abs(f_TF[!msk, 'score']))
tfs <- f_TF %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_TF <- f_TF %>%
  filter(source %in% tfs)

#
f_TF <- f_TF[(abs(f_TF$score) > 0.5),]


c(
  "GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
  "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")

# Plot
ggplot(f_TF, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  # scale_fill_gradient2(low = "yellow3", high = "darkred",
  # scale_fill_gradient2(low = "orange3", high = "yellow3",
  scale_fill_gradient2(low = "yellow1", high = "yellow3",
                       mid = "whitesmoke", midpoint = 0) + 
  theme_bw(base_size = 30) +
  theme(axis.title = element_text(face = "bold", size = 30),
        axis.text.x = 
          element_text(angle = 60, hjust = 1, size =30, face= "bold"),
        axis.text.y = element_text(size =30, face= "bold")
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank()
  ) +
  xlab("Kinase activity")+
  ylab("Activity score") +
  # ggtitle("CLIP1.ROS1 vs GOPC.ROS1")
  # ggtitle("KIF21A.ROS1 vs CLIP1.ROS1")
  ggtitle("CLIP1.ROS1_KD vs CLIP1.ROS1")


##############################################################################################################
# 12 Target from Kinase                                                                                      #
##############################################################################################################

Target_expression <- MOFA_PTM_Ranked_Nodup
colnames(Target_expression) <- c("target","MOFA_score")
MOFA_genes <- Target_expression$target

Kinase_genes <- f_TF$source
net_filter <- omnipath_ptm[omnipath_ptm$tf %in% Kinase_genes,]
net_filter2 <- net_filter[net_filter$target %in% MOFA_genes,]
colnames(net_filter2)


Kinase_Target <- merge(x=net_filter2,y=Target_expression,by='target', all.x = TRUE)
colnames(Kinase_Target) <- c("PTM_target" , "Kinase",  "mor", "likelihood", "MOFA_score")
Kinase_Target2 <- Kinase_Target[,c(2,1,3,5)] %>% arrange(Kinase)

colnames(f_TF) <- c("statistic","Kinase","condition","score","p_value","rnk")
Kinase_Target3 <- merge(x=Kinase_Target2, y=f_TF[,c(2,4)], by='Kinase',all.x = TRUE)
colnames(Kinase_Target3) <- c("Kinase","PTM_target","mor","PTM_score","Kinase_score")
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024/phosphoMS ROS1 and ALK SerThr and Tyr/ROS1")
# write.table(Kinase_Target3,"Targets_PTM_Kinase_GOPC.ROS1 vs CLIP1.ROS1.txt",sep="\t",row.names=FALSE)
# write.table(Kinase_Target3,"Targets_PTM_Kinase_CLIP1.ROS1 vs KIF21A.ROS1.txt",sep="\t",row.names=FALSE)
# write.table(Kinase_Target3,"Targets_PTM_Kinase_CLIP1.ROS1 vs CLIP1.ROS1_KD.txt",sep="\t",row.names=FALSE)



### Load file and filtering PTM by TOP MOFA score
# kinase_Targets <- read.delim("Targets_PTM_Kinase_GOPC.ROS1 vs CLIP1.ROS1.txt")
Kinase_Target3$PTM_score_abs <- abs(Kinase_Target3$PTM_score)

Kinase_Targets_up <- Kinase_Target3[Kinase_Target3$Kinase_score > 3 , ]
Kinase_Targets_up <- Kinase_Targets_up[Kinase_Targets_up$PTM_score > 0 , ]
Kinase_Targets_up$Subtype <- c("up")


Kinase_Targets_down <- Kinase_Target3[Kinase_Target3$Kinase_score < -2 , ]
Kinase_Targets_down <- Kinase_Targets_down[Kinase_Targets_down$PTM_score < 0 , ]
Kinase_Targets_down$Subtype <- c("down")

Kinase_Targets4 <- rbind(Kinase_Targets_down,Kinase_Targets_up)


Kinase_Targets5 <- Kinase_Targets4 %>% 
  group_by(Kinase) %>%
  arrange(desc(PTM_score_abs), .by_group = TRUE)

Kinase_Targets6 <- Kinase_Targets5[,c(1:5,7,6)]
Kinase_Targets7 <- Kinase_Targets6 %>% 
  group_by(Kinase) %>% 
  top_n(n = -6)

Kinase_Targets8 <- Kinase_Targets7 %>% 
  arrange(Kinase_score)

row.names(Kinase_Targets8) <- Kinase_Targets8$PTM
row.names(Kinase_Targets8) <- make.names(Kinase_Targets8$PTM_target,unique=TRUE)

split <- Kinase_Targets8$Kinase



# Circo plots
# col_fun1 = colorRamp2(c(-3, 0, 3), c("yellow3", "white", "darkred"))
# col_fun1 = colorRamp2(c(-3, 0, 3), c("orange3", "white", "yellow3"))
col_fun1 = colorRamp2(c(-3, 0, 3), c("yellow1", "white", "yellow3"))
col_fun2 = colorRamp2(c(-1, 0, 1), c("white", "white", "white"))


### 1st outer layer
mat <- as.matrix(Kinase_Targets8[,"Kinase_score"])
row.names(mat) <- row.names(Kinase_Targets8)
circos.heatmap(mat, split = split, col = col_fun1,
               rownames.side = "outside",
               rownames.col = 1,
               rownames.cex = 1,
               rownames.font = 2)


### 2nd outer layer
# all_facing_options = c("EGFR" , "SRC", "PRKCA", "GSK3B", "MAPK14", "CLK1", "PAK2", "AKT1", "CDK1", "MAPK10", "MAPK3")
# all_facing_options = c("CDK1" , "CDK2", "CSNK2A1", "MAPK1", "EGFR", "SRC")
all_facing_options = c("CDK9" , "MAP3K13", "CDK1", "CDK18", "CDK20", "CDK5","PIK3C3","MAPKAPK5", "PRKACA","RPS6KB1","INSR","EGFR")
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(-9.5, "mm"), 
              all_facing_options[CELL_META$sector.numeric.index],
              facing = "bending.outside", cex = 1.2,
              adj = c(0.5, 0), niceFacing = TRUE)
  
}, bg.border = NA)

### 3rd outer layer
PTM_id <- as.matrix(Kinase_Targets8[,"PTM_score"])
circos.track(ylim = range(PTM_id), panel.fun = function(x, y) {
  y = PTM_id[CELL_META$subset]
  y = y[CELL_META$row_order]
  circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
  circos.points(seq_along(y) - 0.5, y, col = ifelse(y > 0, "red", "blue"), pch=16
  )
}, cell.padding = c(0.02, 0, 0.02, 0))


lgd = Legend(title = "Activity", col_fun = col_fun1)
grid.draw(lgd)


circos.clear()
dev.off()



#Protein selection by highest SD-------------------------------------------------------------------------------------------

#Select proteins with highest SD
# vars <- apply(mat, 1, sd)
# mat_TOP <- mat[vars > quantile(vars, 0.75), ] 
# mat <- as.matrix(mat_TOP)

### Most variable proteins by condition / group
# mat <- prot.10na
mat <- mean.centered.prot.10na

Pval <- function(x)
{
  trans <- aov(x ~ Sample_Plan$Sample_name)
  return(summary(trans)[[1]][1,5])
}

t.mat <- apply(mat, 1, Pval)
t.mat <- t.mat[t.mat <= 1e-25]
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


### Subgroups extraction
hc <- row_order(ht)

for (i in 1:length(hc)){
  if (i == 1) {
    clu <- t(t(row.names(mat[hc[[i]],])))
    proteomic <- cbind(clu, paste("cluster_PHOSPHO_", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(row.names(mat[hc[[i]],])))
    clu <- cbind(clu, paste("cluster_PHOSPHO_", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}


proteomic
colnames(proteomic) <- c("Gene","Cluster_PHOSPHO") 
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper")
# write.table(proteomic,"list_Oncofusion_cluster_PHOSPHO_TOP1019.txt",sep="\t",row.names=FALSE)




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
mat_Condition <- as.matrix(Condition)

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
                Oncofusion = mat_Condition[1,],
                col = list(
                  Oncofusion = c("ALK" = "blue", "ROS1" = "red")
                )
              ))

ht = draw(ht)



#Subgroups extraction------------------------------
hc <- column_order(ht)

for (i in 1:length(hc)){
  if (i == 1) {
    clu <- t(t(colnames(mat[,hc[[i]]])))
    proteomic <- cbind(clu, paste("cluster", i, sep=""))
    colnames(proteomic) <- c("New_name", "Proteomic")
  } else {
    clu <- t(t(colnames(mat[,hc[[i]]])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    proteomic <- rbind(proteomic, clu)
  }
}

proteomic



proteomic1 <- replace(proteomic[,], proteomic[,] == "cluster1", "HLF positive")
proteomic2 <- replace(proteomic1[,], proteomic1[,] == "cluster2", "PBX positive")

rownames(proteomic2) = proteomic2[,"New_name"]

proteomic3  <- subset(proteomic2, select = -c(New_name))
proteomic5B <- t(proteomic3)



#Monti ConsensusCluster-------------------------------------------------------------------------------------------------


rcc = ConsensusClusterPlus(mat,maxK=3,reps=100,pItem=0.8,pFeature=1,title="example2",distance="pearson",
                           innerLinkage="average",finalLinkage ="ward.D",clusterAlg="hc", seed = 123)


?ConsensusClusterPlus
# resICL = calcICL(rcc,title="example")
# rcc3 = ConsensusClusterPlus(mat,maxK=8,reps=100,pItem=0.8,pFeature=1,title="example3",distance="euclidean",clusterAlg="km")

Nclust <- 2
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



# With imputation
colors <- c("red","blue")
column_tree = hclust(as.dist(1-cor(mat_cons2, method="pearson")), method = "ward.D")
column_tree = color_branches(column_tree, k = Nclust, col = colors)
# No imputation
# column_tree = hclust(dist(t(mat_cons2), method = "euclidean"), method = "ward.D")
# column_tree = color_branches(column_tree, k = Nclust, col = colors)


#Cluster for rows
column_tree2 = hclust(as.dist(1-cor(t(mat_cons2), method="pearson")), method = "ward.D")
# column_tree2 = hclust(dist(mat_cons2, method = "euclidean"), method = "ward.D")
column_tree2 = color_branches(column_tree2, k = Nclust)



ht <- Heatmap(mat_cons2, na_col = "white", 
              # col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("white", "lightblue", "blue1", "blue2", "blue3")),
              # col = colorRamp2(c(0.8, 0.85, 0.9, 0.95, 1), c("darkblue", "lightblue", "beige", "red1", "darkred")),
              col = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue", "lightblue", "beige", "red1", "darkred")),
              column_dend_reorder = FALSE,
              column_title = "Consensus Clustering",
              name = "Consensus score", 
              cluster_rows = column_tree2,
              cluster_columns = column_tree,
              column_split = Nclust,
              row_split = Nclust,
              show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 5), 
              column_names_gp = gpar(fontsize = 15), 
              top_annotation = HeatmapAnnotation(Proteomic = anno_block(gp = gpar(fill = colors))),
              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = colors))
                                              
              ) 
              
)


ht = draw(ht)

dev.off()
# Monte Carlo Clustering-------------------------------

mat <- t(prot.10na_imput)

library(M3C)
library(NMF) # loading for aheatmap plotting function
library(gplots) # loading this for nice colour scale
library(ggsci) # more cool colours

pca(mat,legendtextsize = 10,axistextsize = 10,dotsize=2)

res <- M3C(mat, maxK = 10, method = 1, removeplots = FALSE, iters=100, cores=1, seed = 456, clusteralg = "hc")
# res$scores
# res$plots
# M3C

Nclust <- 5
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
Cell_Line <- SamplePlan2["Sample_name",]
mat_CL <- as.matrix(Cell_Line)



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
              column_title = "Unsupervised M3C clustering_ROS1 (k=5)",
              name = "Consensus score", 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              # column_split = 8,
              # show_row_names = FALSE, 
              row_names_gp = gpar(fontsize = 5), 
              column_names_gp = gpar(fontsize = 12), 
              top_annotation = HeatmapAnnotation(
                Oncofusion = mat_G1[1,],
                Cell_Line = mat_CL[1,],
                col = list(
                  Oncofusion = c("ALK" = "blue", "ROS1" = "red"),
                  Cell_Line = c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4",
                                "GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                                "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")
                )
              ))

ht = draw(ht)


proteomic5C <- SamplePlan2[c(4,5),]



#Proteome matrix with proteome classification-----------------------


data_1pept2B <- data_1pept2
sapply(data_1pept2B,class)
data_1pept2B[] <- lapply(data_1pept2B, function(x) as.factor(as.numeric(x)))
sapply(data_1pept2B,class)

proteomic6 <- rbind(proteomic5C,data_1pept2B[, colnames(proteomic5C)])

# proteomic6 <- proteomic6[,-c(4,5)]
# saveRDS(proteomic6, "Full_Surfaceome_matrix_ASP14_CL_sub_info.rds")

# saveRDS(proteomic6, "Full_Surafaceome_matrix_CL_sub_info.rds")



# test1 <- as.data.frame(data_1pept2["INSR","LK1850"])
# test2 <- as.data.frame(proteomic6[c("Proteome_group","INSR"),"LK1850"])
# head(test1)
# head(test2)
# head(proteomic5B[,"LK1850"])






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
names(my_data_PCA2)[names(my_data_PCA2) == 'proteomic10$Sample_name'] <- 'Cell_Line'



g <- ggplot(my_data_PCA2, aes(PC1, PC2, label = c(rownames(my_data_PCA2))
                              , col = Cell_Line, fill = Cell_Line
                              )) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.3) +
  geom_point(size = 5, shape = 21, col = "black") +
  theme_light(base_size=20) +
  
  
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA_ROS1_TOP5000 PTMs")
  

g1 <- g + scale_fill_manual(values = c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                                       "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"))

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
                 Cell_Line = proteomic10$Sample_name)





t <- ggplot(df, aes(x, y, colour = Cell_Line)) +
  geom_point(size = 6,pch = 16) +
  theme_light(base_size=20) +
  
  
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("UMAP_6975 proteins") 


t1 <- t + scale_colour_manual(values = c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                                         "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"))
  
t1



