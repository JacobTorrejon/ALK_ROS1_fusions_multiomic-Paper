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
rna <- readRDS("RNAseq matrix_ROS1_Normalized_filtered_protein coding genes.rds")
# row.names(rna)<-paste(row.names(rna),"RNA",sep="_")

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Proteome")
prot <- readRDS("Proteome matrix_ROS1_Normalized.rds")
# row.names(prot)<-paste(row.names(prot),"PROT",sep="_")

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Phospho/Phospho_2024/phosphoMS ROS1 and ALK SerThr and Tyr/ROS1")
phospho <- readRDS("Phospho matrix_ROS1_Normalized.rds")

# Long_phospho_name <- row.names(phospho)
Short_phospho_name <- sub("\\.","-Phospho:",row.names(phospho))
# Names <- as.data.frame(cbind(Long_phospho_name,Short_phospho_name))
row.names(phospho) <- Short_phospho_name



# Sample Plan
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/Proteome")
# Sample_Plan <- read.delim("Sample_Plan.txt")
Sample_Plan <- read.csv("Sample_Plan2.csv")
row.names(Sample_Plan) <- Sample_Plan$Replica_name
meta <- Sample_Plan[17:34,]
names(meta)[names(meta) == 'Replica_name'] <- 'sample'
names(meta)[names(meta) == 'Sample_name'] <- 'Cell_Line'
rownames(meta) <- meta$sample

meta$Cell_Line <- as.factor(meta$Cell_Line)
meta$Oncofusion <- as.factor(meta$Oncofusion)

#Selecting one group-----------------------------------------------------------------------

# meta_sub <- subset(meta, meta[,"Group"] == "Oligo")
# meta_sub2 <- subset(meta_sub, meta_sub[,"Proteome"] == "Oligo_prot1" | meta_sub[,"Proteome"] == "Oligo_prot2" | meta_sub[,"Proteome"] == "Oligo_prot3")
# 
# 
# list_sub <- row.names(meta_sub2)
# prot_sub <- prot[, colnames(prot) %in% list_sub]
# prot <- prot_sub



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
# rna_p = t(scale(t(rna_p), center = TRUE, scale = FALSE))
# prot_p = t(scale(t(prot_p), center = TRUE, scale = FALSE))
# phospho_p = t(scale(t(phospho_p), center = TRUE, scale = FALSE))


# Centering sample distribution
rna_p = scale(rna_p, center = TRUE, scale = FALSE)
prot_p = scale(prot_p, center = TRUE, scale = FALSE)
phospho_p = scale(phospho_p, center = TRUE, scale = FALSE)

boxplot(rna_p[,1:18])
boxplot(prot_p[,1:18])
boxplot(phospho_p[,1:18])


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


MOFAobject.trained <- run_mofa(MOFAobject, use_basilisk = TRUE)
setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ROS1")
# saveRDS(MOFAobject.trained,"MOFA_RNA_Prot_Phospho_ROS1_ANOVA_pval0.01_5factors.rds")
MOFAobject.trained



#### SECOND PART DOWNSTREAM ANALYSIS----------------------------------------------------------------
### 6. PRE-TRAINED MODEL IMPORTATION

setwd("~/Documents/Curie_Jacob/Proteomic/Oncofusion driven gliomas/MOFA/MOFA_RNA_Prot_NewPhospho/ROS1")
MOFAobject <- readRDS("MOFA_RNA_Prot_Phospho_ROS1_ANOVA_pval0.01_5factors.rds")
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

plot_data_overview(MOFAobject, colors=c("#5BB5B5","#B480D6","#E19B5D"))
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=15)
# grid.table(MOFAobject@cache[["variance_explained"]][["r2_per_factor"]][["group1"]])
var <- as.data.frame(MOFAobject@cache[["variance_explained"]][["r2_per_factor"]][["group1"]])
var2 <- var %>% mutate_if(is.numeric, ~sprintf("%.3f",.)) 
grid.table(var2)
plot_variance_explained(MOFAobject, x="view", y="factor", plot_total = T)



### 8. EXPLORATION OF ASSOCIATION WITH COVARIABLES
## A. Pearson correlation 
correlate_factors_with_covariates(MOFAobject,
                                  # plot = "r",
                                  covariates = c("Cell_Line","Oncofusion"),
)


## B. WITH 1 covariable
# Proteome
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
                  factors = c(4,5), 
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




### 9. Correlations quantification---------------------------------------------------------------------------------------
# A. ANOVA: Factor 1  ~ Group, Group2
summary(aov(MOFAobject@expectations[["Z"]][["group1"]][,1]~ MOFAobject@samples_metadata[["Group2"]]))
summary(aov(MOFAobject@expectations[["Z"]][["group1"]][,1]~ MOFAobject@samples_metadata[["Group"]]))


# ANOVA with Tukey
Anova <- aov(MOFAobject@expectations[["Z"]][["group1"]][,1]~ MOFAobject@samples_metadata[["Group2"]])
summary(Anova)
coefficients(Anova)
Tukey <- TukeyHSD(Anova)
plot(Tukey)
Tukey


## ANOVA Boxplot with significant levels-----------------------------------------
library(ggpubr)

# Group2
df.factor1 <- cbind(df.Group2,MOFAobject@expectations[["Z"]][["group1"]][,1])
colnames(df.factor1) <- c("Group2","value")
df.factor1$Group2 <- factor(df.factor1$Group2, exclude = NULL,
                                                    levels = c("Ewing", "No_Ewing"),
                                                    labels = c("Ewing", "No_Ewing"))


boxplot(value ~ Group2, data = df.factor1,boxwex = 0.5, id.n = Inf,
        lwd = 2, xlab = 'Group2', ylab = 'Factor1', col = c("Ewing" = "pink","No_Ewing" = "navy"), outline=FALSE)

stripchart(value ~ Group2, data = df.factor1, vertical = TRUE,
           method = "jitter", add = TRUE, lwd=6, pch = 20, col = c("Ewing" = "pink","No_Ewing" = "navy"))


my_comparisons <- list( c("Ewing", "No_Ewing"))
ggboxplot(df.factor1, x="Group2", y="value", ylab = "Factor1", xlab = "Group2",
          color = "Group2", palette =c("black","black"),
          # add = "jitter",
          legend = "none",
          fill = c("Ewing" = "pink","No_Ewing" = "navy")) + 
  geom_point(shape=19, size=3) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 




# Group
df.factor1 <- cbind(df.Group,MOFAobject@expectations[["Z"]][["group1"]][,1])
colnames(df.factor1) <- c("Group","value")
df.factor1$Group <- factor(df.factor1$Group, exclude = NULL,
                           levels = c("Ewing", "Immortalized MSC",  "Ovarian cancer","Medulloblastoma",
                                      "Neuroblastoma","Pancreatic cancer","Uveal melanoma","Colon cancer","Mesenchymal stem cells"),
                           labels = c("Ewing", "Immortalized MSC",  "Ovarian cancer","Medulloblastoma",
                                      "Neuroblastoma","Pancreatic cancer","Uveal melanoma","Colon cancer","Mesenchymal stem cells"))


boxplot(value ~ Group, data = df.factor1,boxwex = 0.5, id.n = Inf,
        lwd = 2, xlab = 'Group', ylab = 'Factor1', col = c("Ewing" = "pink1", 
                                                           "Immortalized MSC" = "yellow1",  "Ovarian cancer"= "orange1","Medulloblastoma"="orange3",
                                                           "Neuroblastoma"="greenyellow","Pancreatic cancer"="green2","Uveal melanoma"="green4","Colon cancer"="yellow3","Mesenchymal stem cells"="blue"), outline=FALSE)

stripchart(value ~ Group, data = df.factor1, vertical = TRUE,
           method = "jitter", add = TRUE, lwd=6, pch = 20, col = c("Ewing" = "pink1", 
                                                                   "Immortalized MSC" = "yellow1",  "Ovarian cancer"= "orange1","Medulloblastoma"="orange3",
                                                                   "Neuroblastoma"="greenyellow","Pancreatic cancer"="green2","Uveal melanoma"="green4","Colon cancer"="yellow3","Mesenchymal stem cells"="blue"))


my_comparisons <- list( c("Ewing", "Immortalized MSC"),c("Ewing", "Ovarian cancer"),c("Ewing", "Medulloblastoma"),c("Ewing", "Neuroblastoma"),
                        c("Ewing", "Pancreatic cancer"),c("Ewing", "Uveal melanoma"),c("Ewing", "Colon cancer"),c("Ewing", "Mesenchymal stem cells"))
ggboxplot(df.factor1, x="Group", y="value", ylab = "Factor1", xlab = "Group",
          color = "Group", palette =c("black","black","black","black","black","black","black","black","black"),
          # add = "jitter",
          legend = "none",
          fill = c("Ewing" = "pink1", 
                   "Immortalized MSC" = "yellow1",  "Ovarian cancer"= "orange1","Medulloblastoma"="orange3",
                   "Neuroblastoma"="greenyellow","Pancreatic cancer"="green2","Uveal melanoma"="green4","Colon cancer"="yellow3","Mesenchymal stem cells"="blue"
                   )) + 
  geom_point(shape=19, size=3) + 
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 




# ## B. GLM
# t1 <- data.frame(Rbis = MOFAobject@samples_metadata[["Rbis"]], 
#                  Myc_status = MOFAobject@samples_metadata[["Myc_status"]], 
#                  factor1 = MOFAobject@expectations[["Z"]][["group1"]][,1], 
#                  DKFZ_subT = MOFAobject@samples_metadata[["DKFZ_subtype"]],
#                  DKFZ_subG = MOFAobject@samples_metadata[["DKFZ_subgroup"]])
# 
# t1 <- na.omit(t1)
# 
# # Relapse ~ Factor 1
# t1_mod <- glm(Rbis ~ factor1, 
#               family = binomial(link=logit),
#               data=t1)
# t1_mod0 <- glm(Rbis ~ 1, 
#                family = binomial(link=logit),
#                data=t1)
# 
# anova(t1_mod0,t1_mod, test="Chisq")  
# 
# summary(t1_mod)
# 
# # Relapse ~ DKFZ
# t1_mod <- glm(Rbis ~ DKFZ_subG, 
#               family = binomial(link=logit),
#               data=t1)
# t1_mod0 <- glm(Rbis ~ 1, 
#                family = binomial(link=logit),
#                data=t1)
# 
# anova(t1_mod0,t1_mod, test="Chisq")  
# 
# summary(t1_mod)
# 
# # Relapse ~ Factor1 + DKFZ_subG
# 
# t1_mod <- glm(Rbis ~ factor1 + DKFZ_subG, 
#               family = binomial(link=logit),
#               data=t1)
# t1_mod0 <- glm(Rbis ~ 1, 
#                family = binomial(link=logit),
#                data=t1)
# 
# anova(t1_mod0,t1_mod, test="Chisq")  
# 
# summary(t1_mod)


### 10. TOP WEIGHTS----------------------------------------------------------------------------------
n_factor <- 2

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




# weights$PTM <- row.names(weights)
# lapply(weights,dim)
# 
weights_Prot["NFIB",]
weights_RNA["NFIB",]
weights_Phospho["NFIB",]

weights_RNA["MYC",]



### 11. Samples on top weighted features

n_factor <- 5

plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("RNA expression")) + 
  ggtitle(paste("Top positive weighted RNA expression in factor 1")) +
  scale_color_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                              "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")) +
  scale_fill_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                             "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"))



plot_data_scatter(MOFAobject, 
                  view = "RNA",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("RNA expression")) + 
  ggtitle(paste("Top negative weighted RNA expression in factor 1")) +
  scale_color_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                              "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")) +
  scale_fill_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                             "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"))



plot_data_scatter(MOFAobject, 
                  view = "Proteome",
                  factor = n_factor,  
                  features = 20,
                  sign = "positive",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("Proteome expression")) + 
  ggtitle(paste("Top positive weighted Proteome expression in factor 1")) +
  scale_color_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                              "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")) +
  scale_fill_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                             "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"))





plot_data_scatter(MOFAobject, 
                  view = "Proteome",
                  factor = n_factor,  
                  features = 20,
                  sign = "negative",#positive or negative
                  color_by = "Cell_Line"
) + 
  labs(y= paste("Proteome expression")) + 
  ggtitle(paste("Top negative weighted Proteome expression in factor 1")) +
  scale_color_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                              "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")) +
  scale_fill_manual(values=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                             "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1"))






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

annoCol<-list(Cell_Line=c("GOPC.ROS1"="darkred","GOPC.ROS1_KD"="tomato","CLIP1.ROS1"="yellow3","CLIP1.ROS1_KD"="yellow1",
                          "KIF21A.ROS1"="orange3","KIF21A.ROS1_KD"="orange1")
)

Nclust <- 6

plot_data_heatmap(MOFAobject, 
                  view = "Proteome",
                  factor = 3,
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


## GENERAL


df.class <-  data.frame(IDH1_mutation=df.IDH1_mutation,Grade_Oligo=df.Grade_Oligo,Proteome=df.Proteome)
df.class_order <- df.class %>% rownames_to_column('Sample') %>% arrange(Proteome)

level.order <- c("Oligo_prot1", "Oligo_prot2", "Oligo_prot3" )
df.class_order2 <- df.class_order[order(factor(df.class_order$Proteome, levels = level.order)),]
row.names(df.class_order2) <- df.class_order2$Sample
df.class_order2 <- df.class_order2[,-c(1)]
mat_met <- as.data.frame(t(df.class_order2))

df.Allfactors <- as.data.frame(t(factors.matrix[['group1']][,c(1,2)]))


pheatmap(as.matrix(df.Allfactors),
         # col = colorRampPalette(c("darkblue","cyan", "yellow","orange","red","darkred"))(50),
         annotation_col = df.class_order2,
         cluster_rows = FALSE, cluster_cols = TRUE, 
          cutree_cols = 3,
         show_rownames=TRUE,
         show_colnames=TRUE,
         annotation_colors = annoCol)


### UMAP 
model_umap <- run_umap(
  MOFAobject,
  # factors = "all",
  factors = c(1,2,5),
  # groups = "all",
  n_neighbors = 5,
  min_dist = 0.01)

model_umap <- plot_dimred(MOFAobject, method="UMAP",
                          color_by = "Proteome",
                          # shape_by = "Proteome",
                          dot_size = 3,
                          legend =TRUE)

umap.df <- plot_dimred(MOFAobject, method="UMAP", return_data=TRUE,color_by = "Proteome")


model_umap

umap.df


df <- data.frame(x = umap.df$x,
                 y = umap.df$y,
                 Proteome_groups = umap.df$color_by)


t <- ggplot(df, aes(x, y, colour = Proteome_groups)) +
  geom_point(size = 5,pch = 16) +
  theme_light(base_size=20) +
  
  
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("UMAP_MOFA: RNA + Prot + Phospho") 

t1 <- t + scale_colour_manual(values = c("Oligo_prot1" = "orange1", "Oligo_prot2" = "greenyellow", "Oligo_prot3" = "yellow"))

t1


### TSNE

model_tsne <- run_tsne(
  MOFAobject,
  factors = "all",
  # factors = c(1,2,5),
  perplexity=0.1
  # , theta=0.5, dims=3
  )

model_tsne <- plot_dimred(MOFAobject, method="TSNE"
                          # color_by = "Proteome",
                          # # shape_by = "Proteome",
                          # dot_size = 3,
                          # legend =TRUE
                          )

TSNE.df <- plot_dimred(MOFAobject, method="TSNE", return_data=TRUE,color_by = "Proteome")




TSNE.df


df <- data.frame(x = TSNE.df$x,
                 y = TSNE.df$y,
                 Proteome_groups = umap.df$color_by)


t <- ggplot(df, aes(x, y, colour = Proteome_groups)) +
  geom_point(size = 5,pch = 16) +
  theme_light(base_size=20) +
  
  
  xlab("UMAP1") +
  ylab("UMAP2") +
  ggtitle("UMAP_MOFA: RNA + Prot + Phospho") 

t1 <- t + scale_colour_manual(values = c("Oligo_prot1" = "orange1", "Oligo_prot2" = "greenyellow", "Oligo_prot3" = "yellow"))

t1

##############################
# GSEA
##############################
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(GSA)

setwd("C:/Users/jtorr/Curie_Jacob_03 septiembre/MOFA/GSEA_database")
### Reactome
ReactomeWS <- load(file= "ReactomeWS.RData") 

### Hallmark
Hallmark <- load(file= "Hallmark.RData") 



### OmnipathKinaseSubstrate
# data <-  GSA.read.gmt("omnipath_KinaseSubstrate.gmt")
# data <-  GSA.read.gmt("omnipath_KinaseSubstrate_PTM.gmt")
### Phosphosite Plus
data <-  GSA.read.gmt("PSP_PTM.gmt")
## PTMsigdb from UNIPROT
# data <-  GSA.read.gmt("ptmsigdb_uniprot.gmt") 
# data <-  GSA.read.gmt("Prediction_PTM.gmt") 

gene_names <- unlist(data$genesets, use.names=FALSE)
long_dataframe <- cbind(data$geneset.names,gene_names)
colnames(long_dataframe) <- c("Pathways","Genes")
wide_frame <- as.data.frame(long_dataframe)
head(wide_frame)

wide_frame$Value <- 1
wide_frame2 <- reshape(wide_frame, idvar = "Pathways", timevar = "Genes", direction = "wide")
row.names(wide_frame2) <- wide_frame2$Pathways
row.names(wide_frame2) <- paste("Kinase", row.names(wide_frame2), sep = "_")
# wide_frame3 <-  subset(wide_frame2, select=-c(1) ) # Omnipath
wide_frame3 <-  subset(wide_frame2, select=-c(1,2) ) # PhosphositePlus
colnames(wide_frame3)<-gsub("Value.","",colnames(wide_frame3))
wide_frame3[is.na(wide_frame3)] <- 0

# Pathway_Protein <- as.matrix(mydataframe3)
Pathway_PTM <- as.matrix(wide_frame3)
Pathway_PTM2 <- Pathway_PTM
colnames(Pathway_PTM2) <- gsub(x = colnames(Pathway_PTM2), pattern = "\\:", replacement = ".") 
colnames(Pathway_PTM2) <- gsub(x = colnames(Pathway_PTM2), pattern = "\\-", replacement = ".") 


# Intersect_PTM <- Reduce(intersect, list(colnames(Pathway_PTM2), row.names(phospho)))
# 
# library(VennDiagram)
# PSP_PTM <- colnames(Pathway_PTM2)
# Phospho_MOFA <- row.names(phospho)
# 
# myCol <- c("blue", "red")
# venn.diagram(
#   x = list(Phospho_MOFA, PSP_PTM),
#   category.names = c("Phospho_MOFA", "Omnipath"),
#   filename = 'venn_diagramm_MOFA_DDA vs Omnipath.png',
#   output=TRUE,
# 
#   lwd = 2,
#   lty = 'blank',
#   fill = myCol,
# 
# )


enrichment.parametric <- run_enrichment(MOFAobject,
                                        view = "Proteome", factors = c(1:2),
                                        feature.sets = reactome.matrix,
                                        # feature.sets = hallmark.matrix,
                                        # feature.sets = Pathway,
                                        # feature.sets = Pathway_PTM2,
                                        sign = "negative",
                                        set.statistic = "mean.diff",
                                        min.size = 10,
                                        statistical.test = "parametric",
                                        # statistical.test = "permutation",
                                        # nperm = 200
)




plot_enrichment_heatmap(enrichment.parametric)

plot_enrichment(enrichment.parametric, 
                factor = 1, 
                max.pathways = 15
)




plot_enrichment_detailed(enrichment.parametric, 
                         factor = 1, 
                         max.genes = 8, 
                         max.pathways = 15,
                         text_size = 3
)





Pathways <- cbind(enrichment.parametric[["pval.adj"]],enrichment.parametric[["set.statistics"]])
colnames(Pathways) <- c("Pval_F1","Pval_F2","Pval_F3","Pval_F4","Pval_F5","Pval_F6","Pval_F7","Pval_F8","Pval_F9","Pval_F10",
                        "Score_F1","Score_F2","Score_F3","Score_F4","Score_F5","Score_F6","Score_F7","Score_F8","Score_F9","Score_F10")


Pathways_F1 <- as.data.frame(Pathways[,c("Pval_F1","Score_F1")])
Pathways_F1 <- Pathways_F1[(Pathways_F1$Pval_F1 < 0.001),]
Rank_Pathways_F1 <- Pathways_F1[order(Pathways_F1$Pval_F1, decreasing = FALSE),, drop = FALSE]
# write.table(Rank_Pathways_F1,"Rank_Pathways_F1_all_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F1,"Rank_Pathways_F1_pos_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F1,"Rank_Pathways_F1_neg_weigths_Reactome.txt",sep="\t",row.names=TRUE)

Pathways_F2 <- as.data.frame(Pathways[,c("Pval_F2","Score_F2")])
Pathways_F2 <- Pathways_F2[(Pathways_F2$Pval_F2 < 0.001),]
Rank_Pathways_F2 <- Pathways_F2[order(Pathways_F2$Pval_F2, decreasing = FALSE),, drop = FALSE]
# write.table(Rank_Pathways_F2,"Rank_Pathways_F2_all_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F2,"Rank_Pathways_F2_pos_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F2,"Rank_Pathways_F2_neg_weigths_Reactome.txt",sep="\t",row.names=TRUE)

Pathways_F3 <- as.data.frame(Pathways[,c("Pval_F3","Score_F3")])
Pathways_F3 <- Pathways_F3[(Pathways_F3$Pval_F3 < 0.001),]
Rank_Pathways_F3 <- Pathways_F3[order(Pathways_F3$Pval_F3, decreasing = FALSE),, drop = FALSE]
# write.table(Rank_Pathways_F3,"Rank_Pathways_F3_all_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F3,"Rank_Pathways_F3_pos_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F3,"Rank_Pathways_F3_neg_weigths_Reactome.txt",sep="\t",row.names=TRUE)

Pathways_F6 <- as.data.frame(Pathways[,c("Pval_F6","Score_F6")])
Pathways_F6 <- Pathways_F6[(Pathways_F6$Pval_F6 < 0.001),]
Rank_Pathways_F6 <- Pathways_F6[order(Pathways_F6$Pval_F6, decreasing = FALSE),, drop = FALSE]
# write.table(Rank_Pathways_F6,"Rank_Pathways_F6_all_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F6,"Rank_Pathways_F6_pos_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F6,"Rank_Pathways_F6_neg_weigths_Reactome.txt",sep="\t",row.names=TRUE)


Pathways_F10 <- as.data.frame(Pathways[,c("Pval_F10","Score_F10")])
Pathways_F10 <- Pathways_F10[(Pathways_F10$Pval_F10 < 0.001),]
Rank_Pathways_F10 <- Pathways_F10[order(Pathways_F10$Pval_F10, decreasing = FALSE),, drop = FALSE]
# write.table(Rank_Pathways_F10,"Rank_Pathways_F10_all_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F10,"Rank_Pathways_F10_pos_weigths_Reactome.txt",sep="\t",row.names=TRUE)
# write.table(Rank_Pathways_F10,"Rank_Pathways_F10_neg_weigths_Reactome.txt",sep="\t",row.names=TRUE)

