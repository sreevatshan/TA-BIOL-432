# dataGenerato.R
# Script to generate a hypothetical dataset

# Load necessary library
library(tibble)

# Step 1: Generate a vector of 100 species names
species_names <- c("Species_A", "Species_B", "Species_C", "Species_D", "Species_E")
species_vector <- sample(species_names, 100, replace = TRUE)

# Step 2: Generate a vector of 100 values for Limb_width
# Using a normal distribution with mean = 5 and sd = 1, ensuring positive values
set.seed(123) # For reproducibility
limb_width <- abs(rnorm(100, mean = 5, sd = 1)) # Absolute to avoid negatives

# Step 3: Generate a vector of 100 values for Limb_length
# Using a normal distribution with mean = 10 and sd = 2, ensuring positive values
limb_length <- abs(rnorm(100, mean = 10, sd = 2))

# Step 4: Generate a vector of 100 values for Observer
observer_names <- c("Alice", "Bob", "Charlie")
observer_vector <- sample(observer_names, 100, replace = TRUE)

# Step 5: Combine all vectors into a data.frame or tibble
measurements <- tibble(
  Species = species_vector,
  Limb_width = limb_width,
  Limb_length = limb_length,
  Observer = observer_vector
)

# Step 6: Export the data.frame to a CSV file
write.csv(measurements, "measurements.csv", row.names = FALSE)
