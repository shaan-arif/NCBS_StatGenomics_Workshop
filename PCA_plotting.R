install.packages("ggalluvial")
library(ggalluvial)
library(readxl)
install.packages("factoextra")
library(factoextra)


sheet_a <- read_excel(
  "filepath/DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "a"
)

sheet_c <- read_excel(
  "filepath/DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "c"
)

sheet_d <- read_excel(
  "filepath/DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "d"
)

sheet_e <- read_excel(
  "filepath/DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "e"
)


## Code for Boxplot:

library(readxl)
DNA_Signature_Data_PCAWG_and_HMF <- read_excel("filepath/DNA_Signature_Data_PCAWG_and_HMF.xlsx", 
                                               sheet = "a")
library(dplyr)

DNA_Signature_Data_PCAWG_and_HMF <-
  DNA_Signature_Data_PCAWG_and_HMF %>%
  mutate(dataset = case_when(
    grepl("^HMF", Sample_ID) ~ "Primary Tumour",
    grepl("^SA",  Sample_ID) ~ "Metastatic Tumour",
    TRUE ~ "Other"
  ))
library(ggplot2)

ggplot(DNA_Signature_Data_PCAWG_and_HMF,
       aes(x = Tumor_location, fill = dataset)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Tumor location",
    y = "Number of samples",
    fill = "Dataset",
    title = "Primary Vs Metastatic tumour samples across tumor locations"
  )


breast_sample_ids <- read_excel(
  "D:/Nathan/PhD_Life/NCBS_Workshop/DNA_Signature_Data_PCAWG_and_HMF - Copy.xlsx",
  sheet = "a"
) %>%
  filter(Tumor_location == "Breast") %>%
  pull(Sample_ID)

breast_sample_ids

###PCA

# Install required packages if not already installed
# install.packages(c("readxl", "tidyverse", "FactoMineR", "factoextra"))

library(readxl)
library(tidyverse)
library(FactoMineR)
library(factoextra)

# -----------------------------
# 1. Load the data
# -----------------------------
df <- read_excel("D:\\Nathan\\PhD_Life\\NCBS_Workshop\\new_data.xlsx")

# -----------------------------
# 2. Prepare data
# Remove Sample_ID and keep numeric columns
# -----------------------------
df_pca <- df %>% 
  column_to_rownames(var = "Sample_ID")

# Ensure all columns are numeric
df_pca <- df_pca %>% mutate(across(everything(), as.numeric))

# -----------------------------
# 3. Run PCA
# scale = TRUE standardizes variables
# -----------------------------
res.pca <- PCA(df_pca, scale.unit = TRUE, graph = FALSE)

# -----------------------------
# 4. Visualizations
# -----------------------------

# Scree plot (variance explained)
fviz_eig(res.pca)

# PCA individuals plot
fviz_pca_ind(res.pca,
             geom.ind = "point",
             pointshape = 21,
             fill = "steelblue",
             repel = TRUE)


##############FINALRESULTS##################

##PCA###
pca_df$Breast <- ifelse(pca_df$Sample_ID %in% breast_sample_ids, "Breast", "Other")
pca_plot_df <- pca_df %>%
  select(Sample_ID,
         Breast,
         PC1 (37.4%),
         PC2 (8%)`) %>%
              rename(
                PC1 = PC1 (37.4%),
                PC2 = PC2 (8%)
                )
library(ggplot2)

ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Breast)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA Plot (PC1 vs PC2)",
    x = "PC1 (37.4%)",
    y = "PC2 (8%)",
    color = "Group"
    ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "right"
    
### Features contributing to PC1 & PC2
library(tidyverse)

load_df <- as.data.frame(loadings)

load_df$Feature <- rownames(load_df)

top_PC1 <- load_df %>% 
  arrange(desc(PC1)) %>% 
  slice(1:20)
ggplot(top_PC1, aes(x = reorder(Feature, PC1), y = PC1)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top Features Contributing to PC1",
       y = "PC1 Loading",
       x = "Feature")

ggplot(load_df, aes(x = reorder(Feature, PC2), y = PC2,
                    fill = PC2 > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Large/clustered (negative)",
                               "Small/simple (positive)")) +
  labs(title = "PC2 Loadings Show SV-Type Separation",
       y = "PC2 Loading",
       x = "Feature",
       fill = "Interpretation")
    
#PCA Loadings

plot(pca_res, type = "l")
kmeans_res <- kmeans(pca_res$x[, 1:5], centers = 4)
loadings
    
#Total PCA Components

plot(pca_res, type = "l")
res.pca$eig

library(tidyverse)

load_df <- as.data.frame(loadings)
load_df$Feature <- rownames(load_df)

####Enriched features in PC1 space

top_PC1 <- load_df %>% 
  arrange(desc(PC1)) %>% 
  slice(1:20) 

ggplot(top_PC1, aes(x = reorder(Feature, PC1), y = PC1)) +
  geom_col() +
  coord_flip() +
  labs(title = "Top Features Contributing to PC1",
       y = "PC1 Loading",
       x = "Feature")

####Enriched features in PC2 space

ggplot(load_df, aes(x = reorder(Feature, PC2), y = PC2,
                    fill = PC2 > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Large/clustered (negative)",
                               labs(title = "PC2 Loadings Show SV-Type Separation",
                                    y = "PC2 Loading",
                                    x = "Feature",
                                    fill = "Interpretation")
                               

##Weird PCA (We don't know)
ggplot(pca_df, aes(PC1, PC2)) +
geom_point(size = 2, alpha = 0.7) +
labs(title = "PCA: PC1 = SV Burden, PC2 = SV Type",
x = "PC1 (SV burden)",
y = "PC2 (small/simple → large/clustered)")
                               