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
library(reshape2)



### 2 IMPORT DATA AND SAMPLE PLAN 
# MATRICES
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
# Sample_Plan <- read.delim("Sample_Plan.txt")
Sample_Plan <- read.csv("Sample_Plan2.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
meta <- Sample_Plan[1:16,]
names(meta)[names(meta) == 'Replica_name'] <- 'sample'
names(meta)[names(meta) == 'Sample_name'] <- 'Cell_Line'
rownames(meta) <- meta$sample

meta$Cell_Line <- as.factor(meta$Cell_Line)
meta$Oncofusion <- as.factor(meta$Oncofusion)




### 3. PRE-TREATMENT----------------------------------
#Keep same patients for all matrices
t <- Reduce(intersect, list(colnames(rna), colnames(prot), colnames(phospho)))

rna_p <- rna[, t]
prot_p <- prot[, t]
phospho_p <- phospho[, t]
meta_p <- meta[t,]

#matrix creation
rna_p <- data.matrix(rna_p, rownames.force = NULL )
prot_p <- data.matrix(prot_p, rownames.force = NULL )
phospho_p <- data.matrix(phospho_p, rownames.force = NULL)


# Centering feature distribution
# prot_p = t(scale(t(prot_p), center = TRUE, scale = FALSE))
# phospho_p = t(scale(t(phospho_p), center = TRUE, scale = FALSE))


# Centering sample distribution
rna_p = scale(rna_p, center = TRUE, scale = FALSE)
prot_p = scale(prot_p, center = TRUE, scale = FALSE)
phospho_p = scale(phospho_p, center = TRUE, scale = FALSE)

boxplot(rna_p[,1:16])
boxplot(prot_p[,1:16])
boxplot(phospho_p[,1:16])


### 4. FEATURE SELECTION (in case of large number of features)

# By standard deviation (unsupervised selection)

# std.dev <- function(x)
# {
#   return(sd(x, na.rm=TRUE))
# }
# 
# t.prot <- apply(prot_p, 1, std.dev)
# t.prot <- t.prot[t.prot >= 0.15]
# prot_p <- prot_p[names(t.prot),]
# 
# t.phospho <- apply(phospho_p, 1, std.dev)
# t.phospho <- t.phospho[t.phospho >= 0.35]
# phospho_p <- phospho_p[names(t.phospho),]


# Or by ANOVA test (supervised selection)
Pval <- function(x)
{
  trans <- aov(x ~ meta_p$Cell_Line)
  # trans <- aov(x ~ meta_p$Group)
  return(summary(trans)[[1]][1,5])
}


t.rna <- apply(rna_p, 1, Pval)
t.rna <- t.rna[t.rna <= 0.001]
rna_p <- rna_p[names(t.rna),]

t.prot <- apply(prot_p, 1, Pval)
t.prot <- t.prot[t.prot <= 0.01]
prot_p <- prot_p[names(t.prot),]

t.phospho <- apply(phospho_p, 1, Pval)
t.phospho <- t.phospho[t.phospho <= 0.01]
phospho_p <- phospho_p[names(t.phospho),]



### FIRST PART: TRAINNING A MODEL------------------------------------------------------------
### 5. CREATION MOFA OBJECT
data <- list(RNA = rna_p, Proteome = prot_p, Phospho = phospho_p)
mobj <-  create_mofa(data)

#Data options
data_opts <- get_default_data_options(mobj)
data_opts$scale_views = TRUE

#Model options
model_opts <- get_default_model_options(mobj)
model_opts$num_factors <- 4
#Training options
train_opts <- get_default_training_options(mobj)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$maxiter <- 1000
train_opts$drop_factor_threshold <- 0.0001

#Options in model
MOFAobject <- prepare_mofa(mobj,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)


MOFAobject.trained <- run_mofa(MOFAobject, use_basilisk = TRUE)
# setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho")
# saveRDS(MOFAobject.trained,"MOFA_RNA_Prot_NewPhospho_ALK_ANOVA_pval0.01_2factors.rds")
MOFAobject.trained



#### SECOND PART DOWNSTREAM ANALYSIS----------------------------------------------------------------
### 6. PRE-TRAINED MODEL IMPORTATION
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ALK")
MOFAobject <- readRDS("MOFA_RNA_Prot_NewPhospho_ALK_ANOVA_pval0.01_2factors.rds")
# MOFAobject <- MOFAobject.trained
samples_metadata(MOFAobject) <- meta_p
MOFAobject

#Getting model outputs for further graphical representation

df.Group <- MOFAobject@samples_metadata[,c("sample","Cell_Line")]
rownames(df.Group) <- df.Group$sample; df.Group$sample <- NULL

df.Oncofusion <- MOFAobject@samples_metadata[,c("sample","Oncofusion")]
rownames(df.Oncofusion) <- df.Oncofusion$sample; df.Oncofusion$sample <- NULL

factors.frame <- get_factors(MOFAobject, as.data.frame = T)
factors.matrix <- get_factors(MOFAobject)


### 7. MODEL VERIFICATION GRAPHS----------------------------------------------------------------------------------

plot_data_overview(MOFAobject, colors=c("#5BB5B5","#B480D6","#E19B5D"), show_dimensions = TRUE)
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=75)
# grid.table(MOFAobject@cache[["variance_explained"]][["r2_per_factor"]][["group1"]])
var <- as.data.frame(MOFAobject@cache[["variance_explained"]][["r2_per_factor"]][["group1"]])
var2 <- var %>% mutate_if(is.numeric, ~sprintf("%.3f",.)) 
grid.table(var2)
plot_variance_explained(MOFAobject, x="view", y="factor", plot_total = T)


? plot_data_overview

### 8. EXPLORATION OF ASSOCIATION WITH COVARIABLES
## A. Pearson correlation 
correlate_factors_with_covariates(MOFAobject,
                                  # plot = "r",
                                  covariates = c("Cell_Line","Oncofusion"),
)


## B. WITH 1 covariable
# Proteome
MOFAobject@samples_metadata$Cell_Line <- factor(MOFAobject@samples_metadata$Cell_Line, exclude = NULL,
                                                levels = c("CCDC88A:ALK","CCDC88A:ALK_KD","PPP1CB:ALK","PPP1CB:ALK_KD"),
                                                labels = c("CCDC88A:ALK","CCDC88A:ALK_KD","PPP1CB:ALK","PPP1CB:ALK_KD"))

p <- plot_factor(MOFAobject, 
                 factors = c(1:2),
                 color_by = "Cell_Line",
                 dot_size = 3,   
                 dodge = T,          
                 legend = T,          
                 add_boxplot =  T,     
                 boxplot_alpha = 0.5  
)

p <- p + 
  scale_color_manual(values=c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082")) +
  scale_fill_manual(values=c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082"))

print(p)

### Boxplot with Pvalues
library(ggpubr)
df.factor1 <- cbind(df.Group,MOFAobject@expectations[["Z"]][["group1"]][,1])
colnames(df.factor1) <- c("Cell_Line","value")

p.bx.subtprot.G3B <- ggplot(df.factor1, aes(x=Cell_Line, y=value)) +  
  # geom_violin(aes(fill=Cell_Line), alpha = 1, lwd = 1.2, show.legend = FALSE) +
  geom_boxplot(aes(group=Cell_Line, fill = Cell_Line),
               width=0.4, show.legend = FALSE, alpha = 0.5,
               # position=position_dodge(width=0.35), pch=21, size=1,  stroke = 1.5
               ) + 
  # geom_point(colour = "#4D4D4D", size = 3)+
  geom_point(aes( fill = Cell_Line), size = 3,pch = 21, show.legend = FALSE) +
  scale_color_manual(values=c("#536B18","#C8E485","#147082","#92DFEE")) +
  scale_fill_manual(values=c("#536B18","#C8E485","#147082","#92DFEE")) + 
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2,3.5)+
  geom_signif(comparisons = list(c("PPP1CB:ALK", "PPP1CB:ALK_KD"),c("CCDC88A:ALK", "CCDC88A:ALK_KD")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = -0.55,
              step_increase = 0.6,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 1")) + ggtitle("")


df.factor2 <- cbind(df.Group,MOFAobject@expectations[["Z"]][["group1"]][,2])
colnames(df.factor2) <- c("Cell_Line","value")

p.bx.subtprot.G3B <- ggplot(df.factor2, aes(x=Cell_Line, y=value)) +  
  # geom_violin(aes(fill=Cell_Line), alpha = 1, lwd = 1.2, show.legend = FALSE) +
  geom_boxplot(aes(group=Cell_Line, fill = Cell_Line),
               width=0.4, show.legend = FALSE, alpha = 0.5,
               # position=position_dodge(width=0.35), pch=21, size=1,  stroke = 1.5
  ) + 
  # geom_point(colour = "#4D4D4D", size = 3)+
  geom_point(aes( fill = Cell_Line), size = 3,pch = 21, show.legend = FALSE) +
  scale_color_manual(values=c("#536B18","#C8E485","#147082","#92DFEE")) +
  scale_fill_manual(values=c("#536B18","#C8E485","#147082","#92DFEE")) + 
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2,3.5)+
  geom_signif(comparisons = list(c("PPP1CB:ALK", "PPP1CB:ALK_KD"),c("CCDC88A:ALK", "CCDC88A:ALK_KD")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = -0.0,
              step_increase = 0.1,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 2")) + ggtitle("")


## MOFA_PCA
p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "Cell_Line",
                  # shape_by = "Grade 2021",
                  dot_size = 5,
                  show_missing = T
                  # scale = FALSE
)


p <- p + 
  geom_hline(yintercept=(0.0), linetype="dashed") +
  geom_vline(xintercept=(0.0), linetype="dashed")


p <- p + ggtitle("MOFA_ALK") + 
  # stat_ellipse(geom = "polygon", col = "black", alpha = 0.3) +
  geom_point(size = 5, shape = 21, col = "black") +
  theme_gray(base_size = 25) +
  # theme(base_size = 25,
  #            panel.background = element_rect(fill = "#C2C2C2",colour = "#C2C2C2", size = 0.5, linetype = "solid")
  #                                              ) +
  scale_color_manual(values=c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082")) +
  scale_fill_manual(values=c("CCDC88A:ALK"="#536B18","CCDC88A:ALK_KD"="#C8E485","PPP1CB:ALK_KD"="#92DFEE","PPP1CB:ALK"="#147082"))

print(p)





### 8. TOP WEIGHTS----------------------------------------------------------------------------------
n_factor <- 1

plot_top_weights(MOFAobject,
                 view = "RNA",
                 factor = n_factor,
                 nfeatures = 50,     
                 scale = T,
                 sign = "positive")

plot_top_weights(MOFAobject,
                 view = "Proteome",
                 factor = n_factor,
                 nfeatures = 50,     
                 scale = T,
                 sign = "positive")

plot_top_weights(MOFAobject,
                 view = "Phospho",
                 factor = n_factor,
                 nfeatures = 50,     
                 scale = T,
                 sign = "positive")

plot_top_weights(MOFAobject,
                 view = "RNA",
                 factor = n_factor,
                 nfeatures = 50,     
                 scale = T,
                 sign = "negative")

plot_top_weights(MOFAobject,
                 view = "Proteome",
                 factor = n_factor,
                 nfeatures = 50,     
                 scale = T,
                 sign = "negative")

plot_top_weights(MOFAobject,
                 view = "Phospho",
                 factor = n_factor,
                 nfeatures = 50,     
                 scale = T,
                 sign = "negative")




plot_weights(MOFAobject,
             view = "Proteome",
             factor = 1,
             nfeatures = 50,    # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)



weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 1))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor1, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor1.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 1))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor1, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor1.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 1))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor1, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor1.txt",sep="\t",row.names=TRUE)


weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 2))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor2, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor2.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 2))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor2, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor2.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 2))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor2, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor2.txt",sep="\t",row.names=TRUE)

# weights$PTM <- row.names(weights)
# lapply(weights,dim)
# 
weights_Prot["NFIB",]
weights_RNA["NFIB",]
weights_Phospho["NFIB",]

weights_RNA["MYC",]



### 11. Samples on top weighted features

n_factor <- 1

plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("RNA expression")) + 
  ggtitle(paste("Top positive weighted RNA expression in factor 1")) +
  scale_color_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")) +
  scale_fill_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4"))



plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("RNA expression")) + 
  ggtitle(paste("Top negative weighted RNA expression in factor 1")) +
  scale_color_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")) +
  scale_fill_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4"))



plot_data_scatter(MOFAobject, 
                  view = "Proteome",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("Proteome expression")) + 
  ggtitle(paste("Top positive weighted Proteome expression in factor 1")) +
  scale_color_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")) +
  scale_fill_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4"))




plot_data_scatter(MOFAobject, 
                  view = "Proteome",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("Proteome expression")) + 
  ggtitle(paste("Top negative weighted Proteome expression in factor 1")) +
  scale_color_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")) +
  scale_fill_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4"))





plot_data_scatter(MOFAobject, 
                  view = "Phospho",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("PTM expression")) + 
  ggtitle(paste("Top positive weighted PTM expression in factor 1")) +
  scale_color_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")) +
  scale_fill_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4"))


plot_data_scatter(MOFAobject, 
                  view = "Phospho",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("PTM expression")) + 
  ggtitle(paste("Top negative weighted PTM expression in factor 1")) +
  scale_color_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")) +
  scale_fill_manual(values=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4"))







### 11. MOFA CLASSIFICATION

# One omic with n feactures

annoCol<-list(Cell_Line=c("CCDC88A.ALK"="darkblue","CCDC88A.ALK_KD"="cyan","PPP1CB.ALK_KD"="greenyellow","PPP1CB.ALK"="green4")
)

Nclust <- 3

plot_data_heatmap(MOFAobject, 
                  view = "RNA",
                  factor = 1,
                  features = 80,
                  denoise = TRUE,
                  cutree_cols = Nclust,
                  cutree_rows = Nclust,
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = TRUE,
                  scale = "row", # none, column
                  main= paste("RNA_Factor1_Top80 Proteins"),
                  annotation_samples = data.frame(Cell_Line=df.Group),
                  annotation_colors = annoCol
)

plot_data_heatmap(MOFAobject, 
                  view = "Proteome",
                  factor = 1,
                  features = 80,
                  denoise = TRUE,
                  cutree_cols = Nclust,
                  cutree_rows = Nclust,
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = TRUE,
                  scale = "row", # none, column
                  main= paste("Proteome_Factor1_Top80 Proteins"),
                  annotation_samples = data.frame(Cell_Line=df.Group),
                  annotation_colors = annoCol
)

plot_data_heatmap(MOFAobject, 
                  view = "Phospho",
                  factor = 1,
                  features = 80,
                  denoise = TRUE,
                  cutree_cols = Nclust,
                  cutree_rows = Nclust,
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = TRUE,
                  scale = "row", # none, column
                  main= paste("Phospho_Factor1_Top80 Proteins"),
                  annotation_samples = data.frame(Cell_Line=df.Group),
                  annotation_colors = annoCol
)

dev.off()


