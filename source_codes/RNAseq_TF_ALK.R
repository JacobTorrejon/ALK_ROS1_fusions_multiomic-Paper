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

library(ComplexHeatmap)
library(dendextend)
library(circlize)


##############################################################################################
### 1 IMPORT RNAseq DATA AND SAMPLE PLAN 
##############################################################################################
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper/MOFA/git_Hub/v2")
rna <- readRDS("RNAseq matrix_ALK_Normalized_filtered_protein coding genes.rds")

Sample_Plan <- read.csv("Sample_Plan.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
meta <- Sample_Plan
names(meta)[names(meta) == 'Replica_name'] <- 'sample'
names(meta)[names(meta) == 'Sample_name'] <- 'Cell_Line'
rownames(meta) <- meta$sample

meta$Cell_Line <- as.factor(meta$Cell_Line)
meta$Oncofusion <- as.factor(meta$Oncofusion)

t <- colnames(rna)
meta_p <- meta[t,]




##############################################################################################################
# 2. COMPUTE TF ACTIVITY BASED ON MOFA SCORE                                                                 #
##############################################################################################################

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
MOFA_F1 <- read.delim("Rank_weights_RNA_Factor1.txt")
MOFA_F1$Gene <- row.names(MOFA_F1)
MOFA_F1B <- MOFA_F1 %>% mutate(Gene = (gsub("_RNA", "", Gene)))
row.names(MOFA_F1B) <- MOFA_F1B$Gene

MOFA_F1C <- MOFA_F1B[(abs(MOFA_F1B$Factor1) > 0.25),]


net <- get_collectri(organism='human', split_complexes=FALSE)
# saveRDS(net,"TF_target_interactions_Omnipath.rds")

TF <-  run_ulm(mat=MOFA_F1C[, 'Factor1', drop=FALSE], net=net, .source='source', .target='target',
               .mor='mor', minsize = 10)


### PLOT TOP TF
n_tfs <- 50 # 100

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

# Plot
ggplot(f_TF, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "#C8E485", high = "#536B18", # Supervised
                       mid = "whitesmoke", midpoint = 0) + 
  theme_bw(base_size = 20) +
  theme(axis.title = element_text(face = "bold", size = 30),
        axis.text.x = 
          element_text(angle = 60, hjust = 1, size =30, face= "bold"),
        axis.text.y = element_text(size =30, face= "bold")
        # panel.grid.major = element_blank(), 
        # panel.grid.minor = element_blank()
  ) +
  xlab("Transcription Factor")+
  ylab("Activity score")



##############################################################################################################
# 3. COMPUTE TF ACTIVITY MATRIX                                                                       #
##############################################################################################################



t.rna <- MOFA_F1C$Gene
rna_p <- rna[t.rna,]




# Remove NAs and set row names
counts <- as.matrix(rna_p)
head(counts)


design <- meta_p[,c("sample","Cell_Line")]
colnames(design) <- c("sample","condition")
design


# OmnipathR::collectri()
# net <- get_collectri(organism='human', split_complexes=FALSE)
# net


# Run ulm
sample_acts <- run_ulm(mat=counts, net=net, .source='source', .target='target',
                       .mor='mor',  minsize = 10)

sample_acts



# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  filter(statistic == 'ulm') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()


# Scale per TF
sample_acts_mat <- scale(sample_acts_mat)

# Get top tfs by Pvalue
sample_acts_mat2 <- t(sample_acts_mat)


list_TF <- f_TF$source
sample_acts_mat2 <- sample_acts_mat2[list_TF,]
mat_TF <- as.matrix(sample_acts_mat2)

meta_p2 <- meta_p %>% arrange(Cell_Line)
meta_TF <- as.data.frame(t(meta_p2))
mat_G1 <-  as.matrix(meta_TF["Cell_Line",]) 

### Clustering
colors <- c("black","black","black" ,"black")
Nclust <- 3 #


#Cluster for rows
column_tree2 = hclust(dist(mat_TF, method = "euclidean"), method = "complete")
column_tree2 = color_branches(column_tree2, k = Nclust,col = colors)

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper/TF")
# write.table(mat_TF,"TF_matrix_ALK_TOP50.txt",sep="\t",row.names=TRUE)
# write.table(mat_TF,"TF_matrix_ALK_TOP100.txt",sep="\t",row.names=TRUE)



ht1 <- Heatmap(mat_TF,
               col = colorRamp2(c(-1.5, 0, 1.5), c("#0028A5", "#C2C2C2", "#FF2600")),
               # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               # rect_gp = gpar(col = "white", lwd = 0.1),
               heatmap_legend_param = list(nrow = 1),
               row_title = "TF_TOP100",
               column_title = "Transcription factor activity_ALK",
               name = "Z-score", 
               # column_split = Split,
               # cluster_rows = FALSE,
               cluster_columns = FALSE,
               cluster_rows = column_tree2,
               row_split = Nclust,
               # cluster_columns = column_tree,
               show_row_names = TRUE, 
               show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 8), 
               column_names_gp = gpar(fontsize = 5), 
               top_annotation = HeatmapAnnotation(
                 Cell_Line = mat_G1[1,],
                 col = list(
                   Cell_Line = c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082")
                 )                             
               ))




ht = draw(ht1)




