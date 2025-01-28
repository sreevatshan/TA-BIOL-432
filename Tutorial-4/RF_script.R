##loading the csv file:

library(ggplot2)
library(dplyr)
source("http://bit.ly/theme_pub") # Set custom plotting theme
theme_set(theme_pub())
library(tree)
library(rpart)
library(gbm)
library(randomForest)
library(caret)


can_data <- read.csv("/Users/sreevats/Downloads/Sima-paper_1-16s rRNA-sequences/Cohen_CANCERSEEK_liquid_biopsy_2018_modified.csv", header = T)
str(can_data)

## Identifying the missing summary
missing_summary <- colSums(is.na(can_data))
cat("Missing Data Summary:\n")
print(missing_summary)

##

can_data %>%
  select_if(function(x) any(is.na(x))) %>%
  names()

can_data <- can_data %>%
  mutate (AFP = ifelse(is.na(AFP),0,AFP), Angiopoietin_2 = ifelse(is.na(Angiopoietin_2),0,Angiopoietin_2), AXL = ifelse(is.na(AXL),0,AXL),
          CA_125 = ifelse(is.na(CA_125),0,CA_125), CA_15_3 = ifelse(is.na(CA_15_3),0,CA_15_3), CA19_9 = ifelse(is.na(CA19_9),0,CA19_9),
          CD44 = ifelse(is.na(CD44),0,CD44))

##dimension of data

dim(can_data)

##1804 rows and 42 columns

## identify how many types of tumor is present:

tumor_summary <- table(can_data$Tumor_type)
tumor_summary

##Breast Colorectum  Esophagus      Liver       Lung     Normal      Ovary   Pancreas    Stomach 
##208        388         45         44        104        800         54         93         68 

##seperating the training and test data.

training_data <- can_data[seq(1, nrow(can_data), by = 2), ]
test_data <- can_data[seq(2, nrow(can_data), by = 2), ]

###Training the model with test data

##First converting the "Tumor-level" into factor

training_data <- training_data %>%
  mutate(Tumor_type=as.factor(Tumor_type))

##Construction of tree

train_tree <- tree(Tumor_type ~ ., data = training_data)

plot(train_tree)
text(train_tree, cex=0.7, adj = 0)

summary(train_tree)

##Most influential is IL_8

##using predict to construct a confusion matrix for test data using the training_tree
##before that, making the test_data similar to training data

test_data <- test_data %>%
  mutate(Tumor_type=as.factor(Tumor_type))

##now use predict()

predict_test <- predict(train_tree, test_data, type = "class")

##confusion matrix

cat_data <- data.frame(Obs=test_data$Tumor_type, Pred=predict_test)
table(cat_data)

##Pred
##Obs          Breast Colorectum Esophagus Liver Lung Normal Ovary Pancreas Stomach
##reast         44         25         0     0    9     19     1        3       1
##Colorectum      9        122         0     0   19     25     3        4      15
##Esophagus       4          8         0     0    1      3     0        0       3
##Liver           2          9         0     0    2      4     0        3       4
##Lung            5         13         0     0   17      9     4        1       0
##Normal         13         39         0     0    0    340     1        7       0
##Ovary           6         10         0     0    0      1    10        1       0
##Pancreas        1          7         0     0    3      8     4       25       0
##Stomach         5         10         0     0    8      1     0        2       9

##Finding out the misclassification rate

Correct <- cat_data %>%
  filter(Obs==Pred)
nrow(Correct)/nrow(cat_data)

##0.6286031 - Correct

MisClass<-cat_data %>%
  filter(Obs!=Pred)
nrow(MisClass)/nrow(cat_data)

## 0.3713969 - Misclassification

## PART-3 Randomforest

set.seed(123)

## we have 9 different categories so we need to use atleast 8,

rf_training_data <- randomForest(Tumor_type ~., data = training_data, ntree=100, mtry=7, nodesize=8, importance=TRUE)
rf_training_data


##using predict function on this!

rf_predict <- predict(rf_training_data, test_data)

##confusion matrix

cat_data_rf_training_data <- data.frame(Obs=test_data$Tumor_type, Pred=rf_predict)
table(cat_data_rf_training_data)

#Pred
#Obs          Breast Colorectum Esophagus Liver Lung Normal Ovary Pancreas Stomach
#Breast         80         18         0     0    4      0     0        0       0
#Colorectum      5        183         0     0    0      1     0        0       8
#Esophagus       5         13         0     0    0      0     0        0       1
#Liver           6         12         0     3    1      0     0        0       2
#Lung           11         18         0     0   20      0     0        0       0
#Normal          0          0         0     0    0    400     0        0       0
#Ovary           4          3         0     0    0      1    20        0       0
#Pancreas        0         10         0     0    0      2     0       36       0
#Stomach        10         17         0     0    3      0     0        0       5

##Finding out the misclassification rate for rf_model

Correct <- cat_data_rf_training_data %>%
  filter(Obs==Pred)
nrow(Correct)/nrow(cat_data_rf_training_data)

##0.8281596 - Correct

MisClass<-cat_data_rf_training_data %>%
  filter(Obs!=Pred)
nrow(MisClass)/nrow(cat_data_rf_training_data)

##0.1718404- Misclassification

##plotting the important features

# Extract feature importance
feature_importance <- importance(rf_training_data)

# Convert to a data frame for easier plotting
importance_df <- data.frame(
  Feature = rownames(feature_importance),
  MeanDecreaseAccuracy = feature_importance[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = feature_importance[, "MeanDecreaseGini"]
)

# Sort by Mean Decrease Accuracy
importance_df <- importance_df[order(-importance_df$MeanDecreaseAccuracy), ]

# Plot the feature importance
ggplot(importance_df, aes(x = reorder(Feature, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Feature Importance Based on Random Forest",
    x = "Protein Features",
    y = "Mean Decrease in Accuracy"
  ) +
  theme_minimal()

##sFas and IL8 are most important feature!

#####FInal part

##Create a new column "binary" where it says whether the samples is cancer or normal

binary_data <- can_data %>%
  mutate(
    binary = ifelse(Tumor_type == "Normal", "normal", "cancer")
  )

##now repeating everything!

##separating_test and training data

training_data_binary <- binary_data[seq(1, nrow(binary_data), by = 2), ]
test_data_binary <- binary_data[seq(2, nrow(binary_data), by = 2), ]

##now assigning "binary" column as factor

training_data_binary <- training_data_binary %>%
  mutate(binary=as.factor(binary))

test_data_binary <- test_data_binary %>%
  mutate(binary=as.factor(binary))

##training a rfmodel

rf_model_binary <- randomForest(binary ~., data = training_data_binary, ntree=100, mtry=1, nodesize=8, importance=TRUE)

##using predict function on this!

rf_predict_binary <- predict(rf_model_binary, test_data_binary)

##confusion matrix

cat_data_rf_training_data_binary <- data.frame(Obs=test_data_binary$binary, Pred=rf_predict_binary)
table(cat_data_rf_training_data_binary)

#Pred
#Obs      cancer normal
#cancer    495      7
#normal     12    388

##Finding out the misclassification rate for rf_mode_binary

Correct <- cat_data_rf_training_data_binary %>%
  filter(Obs==Pred)
nrow(Correct)/nrow(cat_data_rf_training_data)

##0.9789357 - Correct

MisClass<-cat_data_rf_training_data_binary %>%
  filter(Obs!=Pred)
nrow(MisClass)/nrow(cat_data_rf_training_data)

##0.0210643 - Misclassification

###Calculating accuracy 

# Create confusion matrix
conf_matrix <- table(cat_data_rf_training_data_binary)

# Calculate accuracy
correct_predictions <- sum(diag(conf_matrix))  # Sum of diagonal elements
total_predictions <- sum(conf_matrix)          # Total number of predictions
accuracy <- correct_predictions / total_predictions

accuracy

### [1] 0.9789357

#### calculating cohen's kappa

library(caret)

# Generate confusion matrix and Cohen's kappa
conf_matrix_kappa <- confusionMatrix(
  data = cat_data_rf_training_data_binary$Pred,  # Predicted labels
  reference = cat_data_rf_training_data_binary$Obs,  # Observed labels
  positive = "cancer"  # Specify the positive class
)

# Extract Cohen's kappa
cohens_kappa <- conf_matrix_kappa$overall["Kappa"]

# Display the results
print(conf_matrix_kappa)  # Full confusion matrix and statistics
cat("Cohen's Kappa:", cohens_kappa, "\n")

#Confusion Matrix and Statistics

#Reference
#Prediction cancer normal
#cancer    495     12
#normal      7    388

#Accuracy : 0.9789          
#95% CI : (0.9673, 0.9873)
#No Information Rate : 0.5565          
#P-Value [Acc > NIR] : <2e-16          

#Kappa : 0.9573     


##plotting the important features

# Extract feature importance
feature_importance_binary <- importance(rf_model_binary)

# Convert to a data frame for easier plotting
importance_df_binary <- data.frame(
  Feature = rownames(feature_importance_binary),
  MeanDecreaseAccuracy = feature_importance_binary[, "MeanDecreaseAccuracy"],
  MeanDecreaseGini = feature_importance_binary[, "MeanDecreaseGini"]
)

# Sort by Mean Decrease Accuracy
importance_df_binary <- importance_df_binary[order(-importance_df_binary$MeanDecreaseAccuracy), ]

# Plot the feature importance
ggplot(importance_df_binary, aes(x = reorder(Feature, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Feature Importance Based on Random Forest",
    x = "Protein Features",
    y = "Mean Decrease in Accuracy"
  ) +
  theme_minimal()

##IL-6 and IL-8

#############

#testing out tuneRF for mtry optimization

# Tune the mtry parameter
tuned_rf <- tuneRF(
  x = training_data[, -which(names(training_data) == "Tumor_type")],  # Predictor variables
  y = training_data$Tumor_type,                                       # Response variable
  stepFactor = 1.5,                                                   # Multiplicative factor for mtry
  improve = 0.01,                                                     # Minimum improvement to continue tuning
  ntreeTry = 500,                                                     # Number of trees for each iteration
  trace = TRUE,                                                       # Display progress
  plot = TRUE                                                         # Plot the error vs. mtry
)

# Print the optimal mtry value
cat("Optimal mtry:", tuned_rf[which.min(tuned_rf[, 2]), 1], "\n")

rf_training_data_optimized <- randomForest(Tumor_type ~., data = training_data, ntree=100, mtry=9, nodesize=8, importance=TRUE)
rf_training_data_optimized

rf_predict_optimized <- predict(rf_training_data_optimized, test_data)

cat_data_rf_training_data_optimized <- data.frame(Obs=test_data$Tumor_type, Pred=rf_predict_optimized)
table(cat_data_rf_training_data_optimized)

##Finding out the misclassification rate for rf_model

Correct__optimized <- cat_data_rf_training_data_optimized %>%
  filter(Obs==Pred)
nrow(Correct__optimized)/nrow(cat_data_rf_training_data_optimized)

##0.8370288 - Correct

MisClass_optimized<-cat_data_rf_training_data_optimized %>%
  filter(Obs!=Pred)
nrow(MisClass_optimized)/nrow(cat_data_rf_training_data_optimized)

##0.1629712- Misclassification


### ~1% increase in accuracy!! 

























