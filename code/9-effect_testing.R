# alpha-diversity ----
{
	
}
# beta-diversity ----
{
	# unbalanced adonis
	{
		unbalanced_community <- master$completewvariance[[5]][, 54:ncol(master$completewvariance[[5]])]
		# run PERMANOVA
		adonis <- list()
		adonis$unbalanced <- adonis2(
			unbalanced_community ~ DTM.km + Catchment + Altitude + Month + Temp + pH + PC1 + PC2 + PC3 + PC4, 
			data         = master$completewvariance[[5]],
			permutations = 999,
			by           = "margin",  # Assess marginal effects (Type III SS)
			method       = "bray"
		)
	}
	## plot
	{
		# 2. Run NMDS (distance = "bray")
		set.seed(123)  # for reproducibility
		nmds_res <- metaMDS(unbalanced_community, distance = "bray", k = 2, trymax = 100)
		
		# 3. Inspect stress to ensure a decent solution
		nmds_res$stress
		
		plot_nmds_by_variable <- function(nmds_obj, data, var_name, main_title = NULL) {
			
			# Extract NMDS coordinates (site scores)
			coords <- nmds_obj$points
			
			# Extract the variable from the data frame
			var_vector <- data[[var_name]]
			
			# Create a blank plot
			plot(coords, type = "n",
				xlab = "NMDS1", ylab = "NMDS2",
				main = ifelse(is.null(main_title), var_name, main_title))
			
			# If numeric, color by gradient; if factor, color by discrete levels
			if (is.numeric(var_vector)) {
				# numeric approach (gradient)
				library(scales)
				color_fun <- col_numeric(
					palette = c(viridis(4)),
					domain  = range(var_vector, na.rm = TRUE)
				)
				point_cols <- color_fun(var_vector)
				
				points(coords, pch = 16, col = point_cols)
				
				# Optional legend with min/max
				legend("topright",
					  legend = c(
					  	paste("Low:", round(min(var_vector, na.rm = TRUE), 2)),
					  	paste("High:", round(max(var_vector, na.rm = TRUE), 2))
					  ),
					  col = c(color_fun(min(var_vector, na.rm = TRUE)),
					  	   color_fun(max(var_vector, na.rm = TRUE))),
					  pch = 16,
					  bty = "n",
					  title = var_name)
				
			} else {
				# factor approach (convex hulls or just color points)
				var_factor <- as.factor(var_vector)
				n_levels   <- nlevels(var_factor)
				
				# color points by factor level
				points(coords, pch = 16, col = as.numeric(var_factor))
				
				# optionally add hulls
				for (i in seq_len(n_levels)) {
					group_points <- coords[var_factor == levels(var_factor)[i], , drop = FALSE]
					if (nrow(group_points) >= 3) {
						hull_idx <- chull(group_points)
						polygon(group_points[hull_idx, ],
							   border = i,
							   lty = 2,
							   col = rgb(
							   	col2rgb(i)[1]/255,
							   	col2rgb(i)[2]/255,
							   	col2rgb(i)[3]/255,
							   	alpha = 0.1
							   ))
					}
				}
				legend("topright",
					  legend = levels(var_factor),
					  col = 1:n_levels,
					  pch = 16,
					  bty = "n",
					  title = var_name)
			}
		}
		
		df_unbalanced <- master$completewvariance[[5]]
		
		par(mfrow = c(2, 2))  # 2x2 grid for plotting
		
		plot_nmds_by_variable(nmds_res, df_unbalanced, "DTM.km", main_title = "DTM.km")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "Catchment", main_title = "Catchment")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "Altitude", main_title = "Altitude")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "Month", main_title = "Month")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "Temp", main_title = "Temp")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "pH", main_title = "pH")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "PC1", main_title = "PC1")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "PC2", main_title = "PC2")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "PC2", main_title = "PC3")
		plot_nmds_by_variable(nmds_res, df_unbalanced, "PC2", main_title = "PC4")
				
		par(mfrow = c(1, 1))
	}
	# nested adonis
	{
		balanced_community <- master$balanced[[5]][, 52:ncol(master$balanced[[5]])]
		ctrl <- how(
			blocks = master$balanced[[5]]$Catchment,
			plots  = Plots(strata = master$balanced[[5]]$Site, type = "free"),
			within = Within(type  = "free"),
			nperm  = 999  # Number of permutations
		)
		adonis$balanced <- adonis2(
			balanced_community ~ DTM.km + Catchment + Altitude + Month + Temp + pH + PC1 + PC2 + PC3 + PC4,
			data         = master$balanced[[5]],
			permutations = ctrl,
			by           = "margin"
		)
		rm(balanced_community)
	}
}
print(adonis$unbalanced)
print(adonis$balanced)
{
	data <- master$completewvariance[[5]]
	data$Month <- as.factor(data$Month)
	plm1 <- glm(
		formula = richness ~ PC1 + PC2 + PC3 + PC4 + Temp + Month + pH,
		family = "poisson" (link = "log"),
		data = data
	)
	bc <- boxcox(plm1)
	lambda <- bc$x[which.max(bc$y)]
	plm2 <- glm(
		formula = ((richness^lambda-1)/lambda) ~ PC1 + PC2 + PC3 + PC4 + Temp + Month + pH,
		family = "poisson" (link = "log"),
		data = data
	)
	nb_mod <- glm.nb(richness ~ Temp + pH + Altitude + Month + PC1 + PC2 + PC3 + PC4, data = data)
	summary(nb_mod)
	gam_mod <- gam(richness ~ s(Temp) + s(pH) + Month + s(PC1) + PC2 + PC3 + PC4, 
				family = poisson(link = "log"), data = data)
	
	summary(gam_mod)
	plot(gam_mod, page = 1)
	# set plot frames
	par(mfrow = c(2,2))
	# draw diagnostic plots
	plot(plm1)
	plot(plm2)
	plot(nb_mod)
	# reset plot frames
	par(mfrow = c(1,1))
}
