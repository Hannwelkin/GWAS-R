# Loading libraries
library(ggplot2)
library(dplyr)
library(MASS)
library(qqman)

# Loading the data
load_data <- function(phenotype_path, genotype_path) {
  phen <- read.csv(phenotype_path, sep = "")
  gen <- read.delim(genotype_path)
  return(list(phen = phen, gen = gen))
}

# Data cleaning
clean_data <- function(phen, gen) {
  # Boxplot for checking outliers in phenotype data
  boxplot.phen <- boxplot(phen$phenos, main = "Boxplot of Phenotype")
  outliers <- boxplot.phen$out
  
  # checking number of rows first
  if (!all(rownames(phen) == rownames(gen))) {
    stop("Phenotype - genotype data are not aligned.")
  }
  
  # Remove outliers from phenotype and genotype data
  phen <- phen[-which(phen$phenos %in% outliers), ]
  gen <- gen[-which(phen$phenos %in% outliers), ]
  
  # Removing missing values in genotype data (assuming '9' is missing data code)
  gen[gen == 9] <- NA
  gen <- na.omit(gen)  # Remove rows with missing genotype data
  
  return(list(phen = phen, gen = gen))
}

# Function to visualize data
visualize_data <- function(phen) {
  # Plotting histogram of phenotype data
  ggplot(phen) +
    aes(x = phenos) +
    geom_histogram(bins = 30L, fill = "#0c4c8a") +
    theme_minimal() +
    labs(x = "BMI/phenotype", y = "Frequency", title = "Phenotype Histogram")
  
  # Boxplot for phenotype outliers
  boxplot.phen <- boxplot(phen$phenos, main = "Boxplot of Phenotype")
  outliers <- boxplot.phen$out
  return(outliers)
}

# Function to process genotypic data
# Remove bad entries (if required)
process_genotypes <- function(gen, bad_rows = NULL, missing_value = 9) {
  if (!is.null(bad_rows)) gen <- gen[-bad_rows, ]
  gen[gen == missing_value] <- NA
  return(gen)
}


# Function to prepare the genotype-phenotype matrix
prepare_gwas_data <- function(gen, phen) {
  phen.numeric <- as.numeric(phen$phenos)
  
  Xa <- data.frame(NA_col = rep(NA, nrow(gen)))
  Xd <- data.frame(NA_col = rep(NA, nrow(gen)))
  
  # Encoding genotypes
  Xa <- ifelse(gen[, seq(1, ncol(gen), 2)] == "A1" & gen[, seq(2, ncol(gen), 2)] == "A1", -1, 
               ifelse(gen[, seq(1, ncol(gen), 2)] == "A1" & gen[, seq(2, ncol(gen), 2)] == "A2", 0, 1))
  Xd <- ifelse(gen[, seq(1, ncol(gen), 2)] == "A1" & gen[, seq(2, ncol(gen), 2)] == "A2", 1, 0)

# Function for GWAS regression and beta calculations
perform_gwas <- function(Xa, Xd, phen) {
  xmatrix <- lapply(1:ncol(Xa), function(i) cbind(1, Xa[, i], Xd[, i]))
  
  beta.list <- lapply(xmatrix, function(m) ginv(t(m) %*% m) %*% t(m) %*% phen)
  
  # Store beta values for each SNP
  beta.mu <- sapply(beta.list, function(x) x[1, ])
  
  return(list(beta.mu = beta.mu, beta.list = beta.list))
}

# Function for F-statistics and p-value calculations
calculate_statistics <- function(phen, beta.list, xmatrix) {
  Y.hat <- sapply(1:length(beta.list), function(i) xmatrix[[i]] %*% beta.list[[i]])
  Y.bar <- mean(phen)
  
  # Calculate Sum of Squares for Model (SSM)
  SSM <- sapply(1:length(beta.list), function(i) sum((Y.hat[, i] - Y.bar)^2))
  
  # Calculate Mean Squares for Model (MSM) and Residuals (SSE)
  MSM <- SSM / 2
  SSE <- sapply(1:length(beta.list), function(i) sum((phen - Y.hat[, i])^2))
  MSE <- SSE / (length(phen) - 3)
  
  # Calculate F-statistics
  F.stat <- MSM / MSE
  
  # Calculate p-values
  Pvalue <- pf(F.stat, df1 = 2, df2 = (length(phen) - 3), lower.tail = FALSE)
  
  return(Pvalue)
}

# Function to visualize GWAS results (Manhattan plot and Q-Q plot)
visualize_gwas_results <- function(Pvalue) {
  # Manhattan plot
  Pvalue.df <- data.frame(Pvalue = Pvalue, chr = gen$chr, Snp = gen$snp, Bp = gen$bp)
  manhattan(Pvalue.df, p = "Pvalue", chr = "chr", bp = "Bp", snp = "Snp", main = "GWAS Manhattan Plot")
  # Q-Q plot
  qq(Pvalue, main = "Q-Q plot of P-values")
}

# Function to identify significant SNPs
identify_significant_snps <- function(Pvalue, alpha = 0.05) {
  significant.snp <- which(Pvalue < alpha)
  return(significant.snp)
}

# Main analysis workflow
run_gwas_analysis <- function(phenotype_path, genotype_path) {
  # Load data
  data <- load_data(phenotype_path, genotype_path)
  phen <- data$phen
  gen <- data$gen
  
  # Clean data
  cleaned_data <- clean_data(phen, gen)
  phen <- cleaned_data$phen
  gen <- cleaned_data$gen
  
  # Visualize data
  outliers <- visualize_data(phen)
  
  # Process genotypes
  gen <- process_genotypes(gen)
  
  # Prepare GWAS data
  gwas_data <- prepare_gwas_data(gen, phen)
  Xa <- gwas_data$Xa
  Xd <- gwas_data$Xd
  phen_numeric <- gwas_data$phen
  
  # Perform GWAS regression and beta calculations
  gwas_results <- perform_gwas(Xa, Xd, phen_numeric)
  beta.mu <- gwas_results$beta.mu
  beta.list <- gwas_results$beta.list
  
  # Calculate F-statistics and p-values
  Pvalue <- calculate_statistics(phen_numeric, beta.list, Xa)
  
  # Visualize GWAS results
  visualize_gwas_results(Pvalue)
  
  # Identify significant SNPs
  significant.snp <- identify_significant_snps(Pvalue)
  
  return(significant.snp)
}

# Example usage (replace with your actual file paths)
phenotype_path <- "Omics99_phenotypes_project1.txt"
genotype_path <- "Omics99_genotypes_project1.txt"

# Run the full analysis
significant_snps <- run_gwas_analysis(phenotype_path, genotype_path)
