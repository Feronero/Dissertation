keyvals$com_factors <- c(
	"Month",
	"pH"
)

# Functions ----
{
	## Data Processing
	process_catchment <- function(catchment, com_factor, tax_level) {
		# 1. Collect data from current catchment and tax. level
		data_subset <- subset(
			master$completewvariance[[tax_level]],
			Catchment == catchment
		)
		community_data <- data_subset[
			,
			54:ncol(data_subset)
		]
		factor_vector <- data_subset[[com_factor]]
		# 3. calculate MDS
		mds_result <- metaMDS(
			comm = community_data,
			distance = "bray",
			k = 2,
			trymax = 100
		)
		# 4. calculate ADONIS
		adonis_result <- adonis2(
			formula = community_data ~ factor_vector
		)
		
		# 5. Output a list of the results of each step
		list(
			data = data_subset,
			community_data = community_data,
			factor_vector = factor_vector,
			MDS = mds_result,
			ADONIS = adonis_result
		)
	}
	### NEW PLOTTING FUNCTION: auto hulls or gradient
	plot_mds_auto <- function(mds_obj, factor_vector, title_text) {
		# set tight margins
		par(
			mar = c(2.5, 2.5, 1.5, 1.5),
			mgp = c(1.5, 0.5, 0)
		)
		
		# Draw blank plot with MDS coordinates
		plot(
			mds_obj$points,
			xlim = range(mds_obj$points[, "MDS1"]),
			ylim = range(mds_obj$points[, "MDS2"]),
			type = "n",
			main = title_text,
			xlab = "MDS1",
			ylab = "MDS2"
		)
		
		# Check if factor_vector is truly numeric or a factor
		if (is.numeric(factor_vector)) {
			# -- CONTINUOUS APPROACH: color gradient
			
			# Create a color function from the 'scales' or 'viridis' package, or define your own
			# Example using 'scales::col_numeric'
			library(scales)  # For col_numeric
			color_fun <- col_numeric(palette = c("blue", "green", "yellow", "red"),
								domain = range(factor_vector, na.rm = TRUE))
			
			point_colors <- color_fun(factor_vector)
			points(mds_obj$points, pch = 16, col = point_colors)
			
			# Optionally add a color legend or colorbar
			# Could do something like a small gradient legend in the margin
			# For brevity, let's just add a simple legend with min & max
			legend("topright",
				  legend = c(
				  	paste("Low:", round(min(factor_vector), 2)),
				  	paste("High:", round(max(factor_vector), 2))
				  ),
				  col = c(color_fun(min(factor_vector)), color_fun(max(factor_vector))),
				  pch = 16,
				  bty = "n",
				  title = "Temp (example)") # Adjust title to suit your variable
			
		} else {
			# -- CATEGORICAL APPROACH: hulls & discrete colors
			n_levels <- nlevels(factor_vector)
			
			# Plot points with distinct colors
			points(
				mds_obj$points,
				pch = 16,
				col = as.numeric(factor_vector)
			)
			
			# Draw convex hulls
			for (i in seq_len(n_levels)) {
				group_points <- mds_obj$points[factor_vector == levels(factor_vector)[i], , drop = FALSE]
				if (nrow(group_points) >= 3) {
					hull_indices <- chull(group_points)
					polygon(
						group_points[hull_indices, ],
						lty = 2,
						border = i,
						col = rgb(
							col2rgb(i)[1]/255,
							col2rgb(i)[2]/255,
							col2rgb(i)[3]/255,
							alpha = 0.1
						)
					)
				}
			}
			
			legend(
				"topright",
				legend = levels(factor_vector),
				col = 1:n_levels,
				pch = 16,
				bty = "n"
			)
		}
		
		# restore default margins
		par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
	}
}
# Execution ----
{
	b_div <- list()
	# use "process_catchment" to fill "b_div"
	for (i in keyvals$catchments) {
		b_div[[i]] <- list()
		for (j in keyvals$com_factors) {
			for (x in seq_len(keyvals$n_levels)) {
				b_div[[i]][[j]][[x]] <- process_catchment(i,	j, x)
			}
		}
	}
	# draw all plots
	for (i in keyvals$catchments) {
		for (j in keyvals$com_factors) {
			for (x in seq_len(keyvals$n_levels)) {
				plot_mds_auto(
					mds_obj = b_div[[i]][[j]][[x]]$MDS,
					factor_vector = b_div[[i]][[j]][[x]]$factor_vector,
					title_text = paste(
						i, j, x,
						"p =", b_div[[i]][[j]][[x]]$ADONIS$`Pr(>F)`[1]
					)
				)
			}
		}
	}
}

# open external plotting window
{
	if (.Platform$OS.type == "windows") {
		windows(width = 10, height = 8)
	} else if (.Platform$OS.type == "unix") {
		quartz(width = 10, height = 8)  # For macOS
	} else {
		x11(width = 10, height = 8)  # For Linux/Unix
	}
}

# adonis2 playtest ----

## subset by catchment
{
	adonis_result <- list()
	for (i in keyvals$catchments) {
		# 1. Collect data from current catchment and tax. level
		data_subset <- subset(
			master[[5]],
			Catchment == i
		)
		community_data <- data_subset[
			,
			54:ncol(data_subset)
		]
		# 3. calculate MDS
		mds_result <- metaMDS(
			comm = community_data,
			distance = "bray",
			k = 2,
			trymax = 100
		)
		# 4. calculate ADONIS
		adonis_result[[i]] <- adonis2(
			formula = community_data ~ Month + Alt_class + Total_Mean.mgL + Temp, 
			data = data_subset,
			method = "bray",
			by = "margin"
		)
	}
	for (i in keyvals$catchments) {
		print(adonis_result[[i]])
	}
}

## all data
{
	# 1. Collect data from current catchment and tax. level
	data_subset <- master[[5]]
	community_data <- master[[5]][ , 54:ncol(data_subset)]
	# 3. calculate MDS
	mds_result <- metaMDS(
		comm = community_data,
		distance = "bray",
		k = 2,
		trymax = 100
	)
	# 4. calculate ADONIS
	adonis_result <- adonis2(
		formula = community_data ~ Lithium_mean.μgL,
		data = data_subset,
		method = "bray",
		by = "margin"
	)
	print(adonis_result)
}
