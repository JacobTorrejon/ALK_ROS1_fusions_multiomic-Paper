##############################################################################################
### 0 LOAD LIBRARIES
##############################################################################################
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
library(reticulate)
use_python("/usr/bin/python3.8")


##############################################################################################
### 1 IMPORT DATA AND SAMPLE PLAN 
##############################################################################################
# MATRICES
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper/MOFA/git_Hub/v2")
rna <- readRDS("RNA_matrix_MOFA_ROS1.rds")
prot <- readRDS("PROT_matrix_MOFA_ROS1.rds")
phospho <- readRDS("PHOSPHO_matrix_MOFA_ROS1.rds")


Sample_Plan <- read.csv("Sample_Plan.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
meta <- Sample_Plan
names(meta)[names(meta) == 'Replica_name'] <- 'sample'
names(meta)[names(meta) == 'Sample_name'] <- 'Cell_Line'
rownames(meta) <- meta$sample

meta$Cell_Line <- as.factor(meta$Cell_Line)
meta$Oncofusion <- as.factor(meta$Oncofusion)

#####################################################################################################
### 2. PUT all MATRICES AND SAMPLE PLAN in the same sample order
#####################################################################################################

t <- Reduce(intersect, list(colnames(rna), colnames(prot), colnames(phospho)))

rna_p <- rna[, t]
prot_p <- prot[, t]
phospho_p <- phospho[, t]
meta_p <- meta[t,]

boxplot(rna_p)
boxplot(prot_p)
boxplot(phospho_p)


################################################################################################
### 3. MOFA TRAINNING A MODEL
################################################################################################
### CREATION MOFA OBJECT
data <- list(RNA = rna_p, Proteome = prot_p, Phospho = phospho_p)
mobj <-  create_mofa(data)

#Data options
data_opts <- get_default_data_options(mobj)
data_opts$scale_views = TRUE

#Model options
model_opts <- get_default_model_options(mobj)
model_opts$num_factors <- 6
#Training options
train_opts <- get_default_training_options(mobj)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$maxiter <- 1000
train_opts$drop_factor_threshold <- 0.00001

#Options in model
MOFAobject <- prepare_mofa(mobj,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)


MOFAobject.trained <- run_mofa(MOFAobject, use_basilisk = FALSE)
# saveRDS(MOFAobject.trained,"MOFA_RNA_Prot_Phospho_ROS1_5factors.rds")
MOFAobject.trained



################################################################################################
### 4. EXPLORATION MOFA FACTORS
################################################################################################
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Figures_paper/MOFA/git_Hub/v2")
MOFAobject <- readRDS("MOFA_RNA_Prot_Phospho_ROS1_5factors.rds")
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

#MODEL VERIFICATION GRAPHS
plot_data_overview(MOFAobject, colors=c("#5BB5B5","#B480D6","#E19B5D"))
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=15)
# grid.table(MOFAobject@cache[["variance_explained"]][["r2_per_factor"]][["group1"]])
var <- as.data.frame(MOFAobject@cache[["variance_explained"]][["r2_per_factor"]][["group1"]])
var2 <- var %>% mutate_if(is.numeric, ~sprintf("%.3f",.)) 
grid.table(var2)
plot_variance_explained(MOFAobject, x="view", y="factor", plot_total = T)




### EXPLORATION OF ASSOCIATION WITH ALK MODEL
MOFAobject@samples_metadata$Cell_Line <- factor(MOFAobject@samples_metadata$Cell_Line, exclude = NULL,
                                                levels = c("GOPC:ROS1","GOPC:ROS1_KD","CLIP1:ROS1","CLIP1:ROS1_KD","KIF21A:ROS1","KIF21A:ROS1_KD"),
                                                labels = c("GOPC:ROS1","GOPC:ROS1_KD","CLIP1:ROS1","CLIP1:ROS1_KD","KIF21A:ROS1","KIF21A:ROS1_KD"))

p <- plot_factor(MOFAobject, 
                 factors = c(1:5),
                 color_by = "Cell_Line",
                 dot_size = 3,   
                 dodge = T,          
                 legend = T,          
                 add_boxplot =  T,     
                 boxplot_alpha = 0.25  
)

p <- p + 
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))

print(p)



### Boxplot with Pvalues
library(ggpubr)

### Factor1
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
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) + 
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2.2,2.0)+
  geom_signif(comparisons = list(c("CLIP1:ROS1", "CLIP1:ROS1_KD")),  
              map_signif_level=TRUE, test = "t.test",
              margin_top = 0.1,
              step_increase = 0.6,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 1")) + ggtitle("")

### Factor2
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
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) + 
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2.2,2.0)+
  geom_signif(comparisons = list(c("CLIP1:ROS1", "GOPC:ROS1"),c("GOPC:ROS1","KIF21A:ROS1")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = -0.18,
              step_increase = 0.20,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 2")) + ggtitle("")


### Factor3
df.factor3 <- cbind(df.Group,MOFAobject@expectations[["Z"]][["group1"]][,3])
colnames(df.factor3) <- c("Cell_Line","value")

p.bx.subtprot.G3B <- ggplot(df.factor3, aes(x=Cell_Line, y=value)) +  
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
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2.2,2.0)+
  geom_signif(comparisons = list(c("CLIP1:ROS1", "CLIP1:ROS1_KD")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = 0.1,
              step_increase = 0.15,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 3")) + ggtitle("")



### Factor4
df.factor4 <- cbind(df.Group,MOFAobject@expectations[["Z"]][["group1"]][,4])
colnames(df.factor4) <- c("Cell_Line","value")

p.bx.subtprot.G3B <- ggplot(df.factor4, aes(x=Cell_Line, y=value)) +  
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
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2.2,2.0)+
  geom_signif(comparisons = list(c("GOPC:ROS1", "GOPC:ROS1_KD")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = 0.2,
              step_increase = 0.15,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 4")) + ggtitle("")### Factor4


### Factor5
df.factor5 <- cbind(df.Group,MOFAobject@expectations[["Z"]][["group1"]][,5])
colnames(df.factor5) <- c("Cell_Line","value")

p.bx.subtprot.G3B <- ggplot(df.factor5, aes(x=Cell_Line, y=value)) +  
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
  geom_hline(yintercept = 0, lty = 2, lwd = 1) + theme_bw() + theme(text = element_text(size = 20, angle = 90)) + ylim(-2.2,2.0)+
  geom_signif(comparisons = list(c("KIF21A:ROS1", "KIF21A:ROS1_KD")),   
              map_signif_level=TRUE, test = "t.test",
              margin_top = 0.1,
              step_increase = 0.15,
              extend_line = 0,
              tip_length = 0.03,
              size = 0.5,
              textsize = 10,
              family = "arial",
              vjust = 0)

print(p.bx.subtprot.G3B + labs(x = "", y = "Factor 5")) + ggtitle("")








## MOFA_PCA
p <- plot_factors(MOFAobject, 
                  factors = c(1,2), # (1,2), (3,4)
                  color_by = "Cell_Line",
                  # shape_by = "Grade 2021",
                  dot_size = 5,
                  show_missing = T
                  # scale = FALSE
)


p <- p + 
  geom_hline(yintercept=(0.0), linetype="dashed") +
  geom_vline(xintercept=(0.0), linetype="dashed")


p <- p + ggtitle("MOFA_ROS1") + 
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.3) +
  geom_point(size = 5, shape = 21, col = "black") +
  theme_gray(base_size = 25) +
  # theme(base_size = 25,
  #       panel.background = element_rect(fill = "#C2C2C2",colour = "#C2C2C2", size = 0.5, linetype = "solid")
  # ) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) 
print(p)




################################################################################################
### 5. MOFA RANKING FEATURES
################################################################################################

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




weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 5))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor5, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor5.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 5))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor5, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor5.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 5))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor5, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor5.txt",sep="\t",row.names=TRUE)




weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 4))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor4, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor4.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 4))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor4, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor4.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 4))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor4, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor4.txt",sep="\t",row.names=TRUE)



weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 3))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor3, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor3.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 3))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor3, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor3.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 3))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor3, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor3.txt",sep="\t",row.names=TRUE)




weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 2))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor2, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor2.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 2))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor2, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor2.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 2))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor2, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor2.txt",sep="\t",row.names=TRUE)



weights_Prot <- as.data.frame(get_weights(MOFAobject, views = "Proteome", factors = 1))
Rank_weights_prot <- weights_Prot[order(weights_Prot$Factor1, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_prot,"Rank_weights_Prot_Factor1.txt",sep="\t",row.names=TRUE)

weights_RNA <- as.data.frame(get_weights(MOFAobject, views = "RNA", factors = 1))
Rank_weights_RNA <- weights_RNA[order(weights_RNA$Factor1, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_RNA,"Rank_weights_RNA_Factor1.txt",sep="\t",row.names=TRUE)

weights_Phospho <- as.data.frame(get_weights(MOFAobject, views = "Phospho", factors = 1))
Rank_weights_Phospho <- weights_Phospho[order(weights_Phospho$Factor1, decreasing = TRUE),, drop = FALSE]
# write.table(Rank_weights_Phospho,"Rank_weights_PTM_Factor1.txt",sep="\t",row.names=TRUE)




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
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))



plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("RNA expression")) + 
  ggtitle(paste("Top negative weighted RNA expression in factor 1")) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))




plot_data_scatter(MOFAobject, 
                  view = "Proteome",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("Proteome expression")) + 
  ggtitle(paste("Top positive weighted Proteome expression in factor 1")) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))






plot_data_scatter(MOFAobject, 
                  view = "Proteome",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("Proteome expression")) + 
  ggtitle(paste("Top negative weighted Proteome expression in factor 1")) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))







plot_data_scatter(MOFAobject, 
                  view = "Phospho",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("PTM expression")) + 
  ggtitle(paste("Top positive weighted PTM expression in factor 1")) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))



plot_data_scatter(MOFAobject, 
                  view = "Phospho",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("PTM expression")) + 
  ggtitle(paste("Top negative weighted PTM expression in factor 1")) +
  scale_color_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                              "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")) +
  scale_fill_manual(values=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                             "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5"))








################################################################################################
### 6. HEATMAP BASED ON TOP MOFA FEATURES
################################################################################################

annoCol<-list(Cell_Line=c("GOPC:ROS1"="#BD3902","GOPC:ROS1_KD"="#FEB799","CLIP1:ROS1"="#8F0A2E","CLIP1:ROS1_KD"="#F78CAA",
                          "KIF21A:ROS1"="#F3AB00","KIF21A:ROS1_KD"="#FFE9B5")
)

Nclust <- 2

plot_data_heatmap(MOFAobject, 
                  view = "Phospho",
                  factor = 1,
                  features = 500,
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
                  view = "RNA",
                  factor = 2,
                  features = 500,
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
                  view = "Proteome",
                  factor = 2,
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


