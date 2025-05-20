# Richness ----
{
	plot(richness ~ PC1, data = complete_metals, col = 1)
	points(richness ~ PC2, data = complete_metals, col = 2)
	glm(
		formula = richness ~ PC1,
		data = complete_metals
	)
	data_subset <- master[[5]] %>%
		filter(!if_any(14:30, is.na))
	lm(
		as.formula(paste("richness ~", paste(keyvals$means, collapse = " + "))),
		data = data_subset
	)
}

## Total Mean Metal
{
	## Overall
	{
		richness.TSMC_guass <- glm(
			richness ~ Total_Mean.mgL,
			family = gaussian(link = "identity"),
			data = master[[5]]
		)
		richness.TSMC_poisson <- glm(
			richness ~ Total_Mean.mgL,
			family = poisson(link = "log"),
			data = master[[5]]
		)
		shapiro.test(residuals(richness.TSMC_guass))
		plot(richness.TSMC_guass)
		dispersiontest(richness.TSMC_poisson)
		
	}
	## By Catchment
	{
		
	}
}
## Individual Metals
{
	## Overall
	{
		richness.metals_gauss <- list()
		for (i in paste0(keyvals$metals, "_mean.μgL")) {
			
		}
	}
}

overall <- lm(master[[5]]$richness ~ master[[5]]$Total_Mean.mgL)
summary(overall)

install.packages("AER")
library(AER)

model2 <- glm(
	master[[5]]$richness ~ master[[5]]$Total_Mean.mgL,
	family = poisson(link = "log")
)
summary(model2)

dispersiontest(model2)
print(dispersion_test)

overdispersion_test <- sum(residuals(model2, type = "pearson")^2) / df.residual(model2)
print(overdispersion_test)

par(mfrow = c(2,2))
plot(model2)
par(mfrow = c(1,1))

plot

summary(model2)

anova(model1, model2)

plot(model)
shapiro.test(resid(model))
model2 <- lm(
	master[[5]]$shannon ~ master[[5]]$Total_Mean.mgL
)
shapiro.test(resid(model2))
plot(model2)

summary(model)
summary(model2)

report(shapiro.test(resid(model)))

# a_div ~ elevation
{
	boxplot(
		master[[5]]$richness ~ master[[5]]$Alt_class,
		ylab = "Richness",
		xlab = "Altitude"
	)
	for (i in keyvals$a_indices[-1]) {
		boxplot(
			master[[5]][[i]] ~ master[[5]]$Alt_class,
			ylab = i,
			xlab = "Altitude"
		)
	}
}
# a-div ~ distance to mouth
{
	plot(
		master[[5]]$richness ~ master[[5]]$DTM.km,
		ylab = "Richness",
		xlab = "Distance to Mouth (km)"
	)
	for (i in keyvals$a_indices[-1]) {
		plot(
			master[[5]][[i]] ~ master[[5]]$DTM.km,
			ylab = i,
			xlab = "Distance to Mouth (km)"
		)
	}
}
# a-div ~ metal conc
{
	for (j in paste0(keyvals$metals, "_mean.μgL")) {
		plot(
			master[[5]]$richness ~ master[[5]][[j]],
			ylab = "Richness",
			xlab = j
		)
	}
	for (i in keyvals$a_indices[-1]) {
		for (j in paste0(keyvals$metals, "_mean.μgL")) {
			plot(
				master[[5]]$shannon ~ master[[5]][[j]],
				ylab = "Shannon",
				xlab = j
			)
		}
	}
}



pca_facto <- PCA(PCA_subset[, 12:28], graph = FALSE)

# Get detailed metal contributions
pca_facto$var$contrib # Contributions to PCs (%)
pca_facto$var$cos2    # Representation quality (0-1)

# Enhanced visualization
fviz_pca_biplot(pca_facto, repel = TRUE, col.var = "red")


# Extract scores and loadings
scores <- as.data.frame(PCA_result$x)
loadings <- as.data.frame(PCA_result$rotation)

# Create biplot
ggplot() +
	# Plot samples
	geom_point(data = scores, aes(x = PC1, y = PC2), alpha = 0.6) +
	# Plot metals (loadings)
	geom_segment(data = loadings,
			   aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5),
			   arrow = arrow(length = unit(0.2, "cm")), 
			   color = "red") +
	geom_text(data = loadings,
			aes(x = PC1*5, y = PC2*5, label = rownames(loadings)),
			color = "red", hjust = 0.5) +
	labs(title = "PCA Biplot", x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
		y = paste0("PC2 (", round(variance_explained[2], 1), "%)"))

	master_df <- master[[5]]
	
		master[[5]][ , c(3:10, 12:28, 46, 47)] %>%
		mutate(across(where(is.numeric), scale))
	comm <- master[[5]][, 54:ncol(master[[5]])]
	common_samples <- intersect(env$SampleID, master_df$SampleID)
	env_complete <- env_complete %>% filter(SampleID %in% common_samples)
	comm_complete <- comm[master_df$SampleID %in% common_samples, ]
	# define env. variable permutation structure
	ctrl <- how(
		blocks = master[[5]]$Catchment,
		plots  = Plots(strata = master[[5]]$Site, type = "free"),
		within = Within(type  = "free"),
		nperm  = 999  # Number of permutations
	)
	# run PERMANOVA
	res <- adonis2(
		comm ~ metalConc, 
		data         = meta_df,
		permutations = ctrl,
		method       = "bray"
	)
}

{
	set.seed(123)
	# select relevant columns
	env <- master[[5]][ , c(3:10, 12:28, 46, 47)] %>%
		mutate(across(where(is.numeric), scale))
	
	# Community data: columns 54 onward
	comm <- master[[5]][, 54:ncol(master[[5]])]
	
	# Drop rows with any NA in Cu, Zn, Mg, etc.
	env_complete <- env %>%
		drop_na(Cu, Zn, Mg)  
	
	# Match community rows accordingly
	common_samples <- intersect(env_complete$SampleID, master_df$SampleID)
	env_complete <- env_complete %>% filter(SampleID %in% common_samples)
	comm_complete <- comm[master_df$SampleID %in% common_samples, ]
}

# Example subset: select relevant columns
env <- master[[5]] %>%
	select(SampleID, River, Site, Month, Cu, Zn, Mg, Temperature, pH)

# Community data: columns 54 onward
comm <- master_df[, 54:ncol(master_df)]

# Drop rows with any NA in Cu, Zn, Mg, etc.
env_complete <- env %>%
	drop_na(Cu, Zn, Mg)  

# Match community rows accordingly
common_samples <- intersect(env_complete$SampleID, master_df$SampleID)
env_complete <- env_complete %>% filter(SampleID %in% common_samples)
comm_complete <- comm[master_df$SampleID %in% common_samples, ]