# Run PCA
pca_result <- princomp(z_scaled_data, cor = FALSE)

# Print PCA summary
summary(pca_result)

# Extract and print eigenvalues
eigenvalues <- pca_result$sdev^2
cat("Eigenvalues:\n", eigenvalues, "\n")

# Explore PCA object structure
names(pca_result)

# Eigenvectors (Loadings)
eigenvectors <- pca_result$loadings
cat("Eigenvectors (Loadings):\n")
print(eigenvectors)

# Principal Component Scores
pc_scores <- pca_result$scores
cat("Principal Component Scores:\n")
head(pc_scores)

# Eigenvalues (variance explained by each PC)
cat("Proportion of Variance Explained:\n")
proportion_variance <- eigenvalues / sum(eigenvalues)
print(proportion_variance)

# Display eigenvector loadings for the first few PCs
print(eigenvectors[, 1:4])


