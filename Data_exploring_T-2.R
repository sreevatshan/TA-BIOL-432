library(tidyverse)
library(ggplot2)
library(moments)
#importing the bird dataset
bird_data <- read.csv("/Users/sreevats/Downloads/Sima-paper_1-16s rRNA-sequences/BirdBehaviour.csv")

##Inspecting the structure of data

str(bird_data)

##Extracting the numeric column for modification

numeric_data <- bird_data %>% select(where(is.numeric))

## Identifying the missing summary
missing_summary <- colSums(is.na(bird_data))
cat("Missing Data Summary:\n")
print(missing_summary)

##Calculating skewness and kurtosis

skewness <- apply(numeric_data, 2, skewness)
kurtosis <- apply(numeric_data, 2, kurtosis)

##Z-scaling 

z_scaled_data <- as.data.frame(scale(numeric_data))
head(z_scaled_data)
summary(z_scaled_data)

# Combine z-scaled data with species
z_scaled_data_with_species <- bird_data %>%
  select(Species) %>%
  bind_cols(z_scaled_data)

# Melt the data for visualization
melted_data <- z_scaled_data_with_species %>%
  pivot_longer(cols = -Species, names_to = "Behavior", values_to = "Z_Score")

# Multi-panel plot
ggplot(melted_data, aes(x = Species, y = Z_Score, fill = Behavior)) +
  geom_boxplot() +
  labs(title = "Behavioral Measurements by Species (Z-Scaled Data)",
       x = "Species", y = "Z-Score") +
  theme_minimal() +
  theme(legend.position = "top")

# Calculate correlation matrix
correlation_matrix <- cor(z_scaled_data, use = "complete.obs")

cor_data <- as.data.frame(as.table(correlation_matrix))

# Plot with ggplot2
ggplot(cor_data, aes(Var1, Var2, fill = Freq)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab", name = "Correlation") +
  geom_text(aes(label = round(Freq, 2)), color = "black", size = 3) +
  labs(title = "Correlation Matrix of Z-Scaled Behavioral Data", 
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))