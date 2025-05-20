# Principal Component Analysis On All Metals ----
{
	## Unbalanced PCA on all paired EA sites with complete metal records
	{
		PCA <- list()
		set.seed(123)
		# define varying_metal_cols
		varying_metal_cols <- master$completewvariance[[1]][, names(master$completewvariance[[1]]) %in% keyvals$medians]
		# run PCA
		PCA$unbalanced <- prcomp(
			varying_metal_cols, 
			center = TRUE,
			scale. = TRUE
		)
		# percentage contributions of each metal to each PC (loading vals squared * 100)
		PCA$unbalanced$percent_contrib <- as.data.frame(
			(PCA$unbalanced$rotation^2) * 100
		)
		colnames(PCA$unbalanced$percent_contrib) <- paste0("PC", 1:ncol(PCA$unbalanced$percent_contrib))
		# variance-weighted contributions
		eigenvalues <- PCA$unbalanced$sdev^2
		PCA$unbalanced$weighted_variance <- as.data.frame(
			(PCA$unbalanced$rotation^2 %*% diag(eigenvalues)) / sum(eigenvalues) * 100
		)
		colnames(PCA$unbalanced$weighted_variance) <- paste0("PC", 1:ncol(PCA$unbalanced$weighted_variance))
		rm(eigenvalues)
		# append balanced PCA scores to each taxonomic level of master$completewvariance
		for (i in 1:keyvals$n_levels) {
			for (j in 1:4) {
				master$completewvariance[[i]][[paste0("PC", j)]] <- PCA$unbalanced$x[, j]
			}
			# Reorder columns to move PCs to the left
			master$completewvariance[[i]] <- master$completewvariance[[i]] %>%
				select(PC1, PC2, PC3, PC4, everything())
		}
	}
	## Balanced PCA on 3 sites per river, - Arsenic, - records from October
	## Heavy filtering required to achieve balance. Worth it?
	{
		balanced_metals <- master$balanced[[1]][, names(master$balanced[[1]]) %in% keyvals$medians]
		PCA$balanced <- prcomp(
			balanced_metals, 
			center = TRUE,
			scale. = TRUE
		)
		variance_explained <- PCA$balanced$sdev^2 / sum(PCA$balanced$sdev^2) * 100
		# percentage contributions of each metal to each PC (loading vals squared * 100)
		PCA$balanced$percent_contrib <- as.data.frame(
			(PCA$balanced$rotation^2) * 100
		)
		colnames(PCA$balanced$percent_contrib) <- paste0("PC", 1:ncol(PCA$balanced$percent_contrib))
		# variance-weighted contributions
		eigenvalues <- PCA$balanced$sdev^2
		PCA$balanced$weighted_variance <- as.data.frame(
			(PCA$balanced$rotation^2 %*% diag(eigenvalues)) / sum(eigenvalues) * 100
		)
		colnames(PCA$balanced$weighted_variance) <- paste0("PC", 1:ncol(PCA$balanced$weighted_variance))
		rm(eigenvalues)
		# append balanced PCA scores to each taxonomic level of master$balanced
		for (i in 1:keyvals$n_levels) {
			for (j in 1:4) {
				master$balanced[[i]][[paste0("PC", j)]] <- PCA$balanced$x[, j]
			}
			# Reorder columns to move PCs to the left
			master$balanced[[i]] <- master$balanced[[i]] %>%
				select(PC1, PC2, PC3, PC4, everything())
		}
	}
}
# Principal Component Analysis on Only Heavy Metals
{
	
}
