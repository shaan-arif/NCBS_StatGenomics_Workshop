library(readxl)
library(dplyr)

signatures <- read_excel(
  "DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "d",
  na = c("", "NA")
)

tumor_meta <- read_excel(
  "DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "a",
  na = c("", "NA")
)

signature_map <- read_excel(
  "DNA_Signature_Data_PCAWG_and_HMF.xlsx",
  sheet = "e"
)

# Add metastatic/primary label
signatures$Label <- ifelse(
  grepl("^HMF", signatures$Donor_ID),
  "Metastatic",
  "Primary"
)
signatures$Label <- factor(signatures$Label)

num_data <- signatures %>% select(where(is.numeric))

num_data <- num_data %>% mutate(across(everything(), ~ as.numeric(.)))
num_data_zero <- num_data %>% replace(is.na(.), 0)


library(dplyr)

model_df <- num_data_zero %>%
  mutate(Label = signatures$Label)

model_df$Label <- factor(model_df$Label)
set.seed(123)
train_idx <- sample(nrow(model_df), 0.8 * nrow(model_df))

train_df <- model_df[train_idx, ]
test_df  <- model_df[-train_idx, ]
library(randomForest)

rf_model <- randomForest(Label ~ ., 
                         data = train_df,
                         importance = TRUE,
                         ntree = 500,
                         proximity = TRUE)

rf_model
pred <- predict(rf_model, newdata = test_df)

conf_matrix <- table(Predicted = pred, Actual = test_df$Label)
conf_matrix
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
accuracy
varImpPlot(rf_model, n.var = 20)

# Check column names
colnames(train_df)
table(signatures$Label, grepl("^HMF", signatures$Donor_ID))

nrow(intersect(train_df, test_df))
table(train_df$Label)
table(test_df$Label)

###Run with stricter method
library(caret)

ctrl <- trainControl(method = "cv", number = 10)

rf_cv <- train(
  Label ~ .,
  data = train_df,
  method = "rf",
  trControl = ctrl
)
## check for cross-validation accuracy
rf_cv

set.seed(99)
train_shuffled <- train_df
train_shuffled$Label <- sample(train_shuffled$Label)

#Permutation test to check for leakage
rf_perm <- randomForest(Label ~ ., data = train_shuffled, ntree = 500)
rf_perm$confusion
rf_perm
#model completely failed on shuffled labels
#so the CV accuacy was correct not due to data leakage

## Run model on Breast cancer samples alone
data_merged <- signatures %>%
  left_join(tumor_meta, by = c("Donor_ID", "Sample_ID"))

#check the data
table(data_merged$Tumor_location)
#filter only breast cancer
breast_data <- data_merged %>% 
  filter(Tumor_location == "Breast")
nrow(breast_data)
head(breast_data)
num_data_b <- breast_data %>% select(where(is.numeric))

num_data_b <- num_data_b %>% mutate(across(everything(), as.numeric))

num_data_b_zero <- num_data_b %>% replace(is.na(.), 0)

#Add labels
breast_df <- cbind(num_data_b_zero, Label = breast_data$Label)
breast_df$Label <- factor(breast_df$Label)


#Create test train split
set.seed(123)

train_idx <- sample(seq_len(nrow(breast_df)), size = 0.7 * nrow(breast_df))

train_breast <- breast_df[train_idx, ]
test_breast  <- breast_df[-train_idx, ]


#train Random forest
library(randomForest)

rf_breast <- randomForest(
  Label ~ .,
  data = train_breast,
  ntree = 2000,
  importance = TRUE
)

rf_breast


#predict on test set
pred_breast <- predict(rf_breast, newdata = test_breast)

conf_breast <- table(Predicted = pred_breast, Actual = test_breast$Label)
conf_breast

#Accuracy on test set
accuracy_breast <- sum(diag(conf_breast)) / sum(conf_breast)
accuracy_breast

#Extract top contributing signatures
importance_breast <- importance(rf_breast)
varImpPlot(rf_breast, n.var = 20, main = "Top Signatures Driving Breast Cancer- Primary & Metastatic Classification")


####################################################
#Logistic Regression with Caret
library(caret)

# Set up 10-fold CV
ctrl <- trainControl(method = "cv", number = 10)

# Fit logistic regression (binomial)
logit_cv <- train(
  Label ~ .,
  data = train_df,
  method = "glm",
  family = binomial,
  trControl = ctrl
)

logit_cv

pred_logit <- predict(logit_cv, newdata = train_df)
confusionMatrix(pred_logit, train_df$Label)

library(pROC)
prob_logit <- predict(logit_cv, newdata = train_df, type = "prob")
roc_obj <- roc(train_df$Label, prob_logit$Metastatic)
auc(roc_obj)
##############################################################

#PCA-Logistic
preproc <- preProcess(train_df[, -1], method = c("center", "scale", "pca"))
train_pca <- predict(preproc, train_df)

logit_pca <- train(
  Label ~ .,
  data = data.frame(Label = train_df$Label, train_pca),
  method = "glm",
  family = binomial,
  trControl = ctrl
)
logit_pca
############################
# Regularized Logistic Regression

logit_lasso <- train(
  Label ~ ., 
  data = train_df,
  method = "glmnet",
  trControl = ctrl,
  tuneLength = 10
)
logit_lasso
##############################

##  PLOTS

library(ggplot2)
library(ggthemes)
library(pheatmap)
library(reshape2)
library(pROC)
library(caret)
library(RColorBrewer)
library(tidyr)


#Top 20 signatures from RF
importance_df <- data.frame(
  Signature = rownames(importance(rf_model)),
  MeanDecreaseGini = importance(rf_model)[, "MeanDecreaseGini"]
) %>% arrange(desc(MeanDecreaseGini)) %>% head(20)

ggplot(importance_df, aes(x = reorder(Signature, MeanDecreaseGini),
                          y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Top 20 Important Signatures (Random Forest)",
       x = "Signature", y = "Mean Decrease in Gini")

# top 20 in breast cancer
importance_b <- data.frame(
  Signature = rownames(importance(rf_breast)),
  MeanDecreaseGini = importance(rf_breast)[, "MeanDecreaseGini"]
) %>% arrange(desc(MeanDecreaseGini)) %>% head(20)

ggplot(importance_b, aes(x = reorder(Signature, MeanDecreaseGini),
                         y = MeanDecreaseGini)) +
  geom_col(fill = "darkred") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Top 20 Important Signatures (Breast Cancer Only)",
       x = "Signature", y = "Mean Decrease in Gini")


#confusion matrix
cm_df <- as.data.frame(conf_matrix)

ggplot(cm_df, aes(Actual, Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 6) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal(base_size = 14) +
  labs(title = "Confusion Matrix (Full Dataset Random Forest)")


cm_b_df <- as.data.frame(conf_breast)

ggplot(cm_b_df, aes(Actual, Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 6) +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_minimal(base_size = 14) +
  labs(title = "Confusion Matrix (Breast Cancer RF)")


#### PCA

pca <- prcomp(num_data_zero, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  Label = signatures$Label
)

var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
ggplot(pca_df, aes(PC1, PC2, color = Label)) +
  geom_point(alpha = 0.7, size = 3) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Primary" = "darkblue", "Metastatic" = "red")) +
  labs(title = "PCA of All Mutational Signatures",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)"))


#heatmap of top 10 signatures
top10 <- importance_df$Signature[1:10]
importance_df
mat <- num_data_zero[, top10]
rownames(mat) <- signatures$Sample_ID

annotation <- data.frame(Label = signatures$Label)
rownames(annotation) <- signatures$Sample_ID


pheatmap(mat,
         annotation_row = annotation,
         color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
         show_rownames = FALSE,
         main = "Top 25 Signatures Heatmap")


# ROC Curve
test_numeric <- ifelse(test_df$Label == "Metastatic", 1, 0)

prob_rf <- predict(rf_model, newdata = test_df, type = "prob")[,2]

roc_obj <- roc(test_numeric, prob_rf)

plot(roc_obj, col = "darkblue", lwd = 4, main = "ROC Curve – RF Model")
auc(roc_obj)

top5 <- importance_df$Signature[1:5]

long_df <- signatures %>% 
  select(all_of(top5), Label) %>% 
  pivot_longer(-Label, names_to = "Signature", values_to = "Value")

ggplot(long_df, aes(Signature, Value, fill = Label)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = c("Primary" = "blue", "Metastatic" = "red")) +
  theme_minimal(base_size = 14) +
  labs(title = "Top 5 Important Signatures: Value Distribution")

### MDS plot
mds <- cmdscale(1 - rf_model$proximity, k = 2)

mds_df <- data.frame(
  Dim1 = mds[,1],
  Dim2 = mds[,2],
  Label = train_df$Label
)

ggplot(mds_df, aes(Dim1, Dim2, color = Label)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Primary" = "blue", "Metastatic" = "red")) +
  labs(title = "Random Forest Proximity-Based MDS Plot")




################## MODELS to Predict Cancer Type #############

# Add cancer type by merging metadata
data_merged <- signatures %>%
  left_join(tumor_meta, by = c("Donor_ID", "Sample_ID"))


#Define multi-calss label
data_merged$Cancer_Type <- factor(data_merged$Tumor_location)

# Convert signatures to numeric (logical → numeric)
num_data <- data_merged %>%
  mutate(across(where(is.logical), ~ as.numeric(.))) %>% 
  mutate(across(where(is.numeric), ~ replace_na(., 0)))


model_df <- num_data %>% 
  filter(!is.na(Cancer_Type))


library(caret)

set.seed(123)
train_idx <- createDataPartition(model_df$Cancer_Type, p = 0.8, list = FALSE)

train_df <- model_df[train_idx, ]
test_df  <- model_df[-train_idx, ]

library(randomForest)

set.seed(123)
rf_cancer <- randomForest(
  Cancer_Type ~ .,
  data = train_df,
  ntree = 500,
  importance = TRUE
)
rf_cancer


pred_test <- predict(rf_cancer, newdata = test_df)

test_results <- confusionMatrix(pred_test, test_df$Cancer_Type)
test_results

ctrl <- trainControl(method = "cv", number = 10)

set.seed(99)
rf_cv <- train(
  Cancer_Type ~ .,
  data = train_df,
  method = "rf",
  trControl = ctrl
)

rf_cv

varImpPlot(rf_cancer, n.var = 20,
           main = "Top Signatures Influencing Cancer-Type Classification")
