setwd("F:/TCGA_all_cancer_genome_analysis/Transcriptome_Profiling/TPM/analysis_no_outliers/protein_coding_dataset/central_moments_analysis")
# Read the gene expression dataset
#data <- read.csv("TPM_proteincoding_norm_data.csv", header = TRUE, row.names = 1)
normaldata <- TPM_proteincoding_norm_data
# Create an empty results list
normal_results_list <- list()

# Iterate through each gene
for (i in 1:nrow(normaldata)) {
  # Extract expression data for the current gene
  normalgene_data <- normaldata[i, ]
  
  # Handle missing values
  normalgene_data <- na.omit(as.numeric(normalgene_data))
  
  # Calculate central moments
  normalgene_mean <- mean(normalgene_data)
  normalgene_variance <- var(normalgene_data)
  normalgene_skewness <- sum((normalgene_data - normalgene_mean)^3) / (length(normalgene_data) * sd(normalgene_data)^3)
  normalgene_kurtosis <- sum((normalgene_data - normalgene_mean)^4) / (length(normalgene_data) * sd(normalgene_data)^4)
  
  # Store results in a data frame
  normal_gene_results <- data.frame(Gene = rownames(normaldata)[i], Mean = normalgene_mean, Variance = normalgene_variance,
                                    Skewness = normalgene_skewness, Kurtosis = normalgene_kurtosis)
  
  # Add results to the results list
  normal_results_list[[i]] <- normal_gene_results
}

# Combine results into a data frame
normal_results <- do.call(rbind, normal_results_list)
write.csv(normal_results, file="normal_results.csv", quote = FALSE)

# Cancer data
cancerdata <- TPM_proteincoding_cancer_data
# Create an empty results list
cancer_results_list <- list()

# Iterate through each gene
for (i in 1:nrow(cancerdata)) {
  # Extract expression data for the current gene
  cancergene_data <- cancerdata[i, ]
  
  # Handle missing values
  cancergene_data <- na.omit(as.numeric(cancergene_data))
  
  # Calculate central moments
  cancergene_mean <- mean(cancergene_data)
  cancergene_variance <- var(cancergene_data)
  cancergene_skewness <- sum((cancergene_data - cancergene_mean)^3) / (length(cancergene_data) * sd(cancergene_data)^3)
  cancergene_kurtosis <- sum((cancergene_data - cancergene_mean)^4) / (length(cancergene_data) * sd(cancergene_data)^4)
  
  # Store results in a data frame
  cancer_gene_results <- data.frame(Gene = rownames(cancerdata)[i], Mean = cancergene_mean, Variance = cancergene_variance,
                                    Skewness = cancergene_skewness, Kurtosis = cancergene_kurtosis)
  
  # Add results to the results list
  cancer_results_list[[i]] <- cancer_gene_results
}

# Combine results into a data frame
cancer_results <- do.call(rbind, cancer_results_list)
write.csv(cancer_results, file="cancer_results.csv", quote = FALSE)

# Calculate skewness ratio
data <- data.frame(normal_results$Skewness, cancer_results$Skewness)

Skewness_result <- apply(data, 1, function(row) {
  max_value <- max(row)
  min_value <- min(row)
  (max_value - min_value) / max_value
})
write.csv(Skewness_result, file="Skewness_result.csv", quote = FALSE)

# Print results
print(result)
