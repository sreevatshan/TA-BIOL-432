# Load necessary library
library(dplyr)

# Read input CSV file (measurements.csv)
data <- read.csv("measurements.csv")

# Check if the required columns are present
required_columns <- c("Species", "Limb_width", "Limb_length", "Observer")
if (!all(required_columns %in% colnames(data))) {
  stop("The input file must contain the columns: Species, Limb_width, Limb_length, Observer")
}

# Define an equation to estimate limb volume
# Using a cylindrical approximation: Volume = π * (radius^2) * height
# Assume radius = Limb_width / 2 and height = Limb_length
data <- data %>%
  mutate(
    Volume = pi * (Limb_width / 2)^2 * Limb_length
  )

# Overwrite the original file with the new data including the Volume column
write.csv(data, "measurements.csv", row.names = FALSE)