# Raw Metal Values ----
{
	# export histograms of raw metal concentrations
	{
		for (i in 1:17) {
			png(
				filename = paste0(
					"outputs/plots/rawmetal_hists/rawmetal_hist_",
					colnames(metal_data$wide[3 + i]),
					".png"
				),
				width = 1080,
				height = 180,
				pointsize = 16,
				res = 100
			)
			par(
				mar = c(0.2, 0, 0.2, 0),
				xaxs = "i",
				yaxs = "i"
			)
			hist(
				metal_data$wide[[3 + i]],
				main = NULL,
				xlab = NULL,
				ylab = NULL,
				axes = FALSE,           # Remove axis framework
				border = "white",       # Optional: remove bar borders
				col = "grey30",         # Optional: set fill color
				panel.first = NULL,     # Remove any default grid
				ann = FALSE,            # Remove all annotations
				xaxt = "n",             # Explicitly remove x-axis
				yaxt = "n",             # Explicitly remove y-axis
				bty = "n"               # Remove bounding box
			)
			box(bty = "n", col = "white")  # Final clean bounding box removal
			dev.off()
		}
		# reset margins
		par(
			mar = c(5.1, 4.1, 4.1, 2.1)
		)
	}
	# Shapiro-Wilk test on all raw metal values
	{
		norm_tests <- list()
		norm_tests$rawmetal <- list()
		for (i in 1:17) {
			x <- metal_data$wide[[3 + i]]
			# if all concentrations identical
			if (length(unique(x)) == 1) {
				norm_tests$rawmetal[[i]] <- data.frame(
					Metal = colnames(metal_data$wide)[3 + i],
					W = NA,
					p.value = NA,
					n = sum(!is.na(x)),
					Note = "Constant values"
				)
				# if concentrations show variance
			} else {
				test <- shapiro.test(x)
				norm_tests$rawmetal[[i]] <- data.frame(
					Metal = colnames(metal_data$wide)[3 + i],
					W = test$statistic,
					p.value = test$p.value,
					n = sum(!is.na(x)),
					Note = ""
				)
			}
		}
		rm(test)
		norm_tests$rawmetal <- do.call(rbind, norm_tests$rawmetal)
		# export table
		write_xlsx(
			x = norm_tests$rawmetal,
			path = "outputs/tables/distest_rawmetal.xlsx"
		)
	}
}

# Median Metal Values ----
{
	metal_points$complete <- metal_points$paired %>%
		filter(if_all(14:30, ~ !is.na(.x)))
	
	## export histograms of median metal concentrations
	{
		for (i in 1:17) {
			png(
				filename  = paste0(
					"outputs/plots/medianmetal_hists/medianmetal_hist_",
					colnames(metal_points$complete[13 + i]),
					".png"
				),
				width     = 1080,
				height    = 180,
				pointsize = 16,
				res       = 100
			)
			par(
				mar  = c(0.2, 0, 0.2, 0),
				xaxs = "i",
				yaxs = "i"
			)
			hist(
				metal_points$complete[[13 + i]],
				main = NULL,
				xlab = NULL,
				ylab = NULL,
				axes = FALSE,           # Remove axis framework
				border = "white",       # Optional: remove bar borders
				col = "grey30",         # Optional: set fill color
				panel.first = NULL,     # Remove any default grid
				ann = FALSE,            # Remove all annotations
				xaxt = "n",             # Explicitly remove x-axis
				yaxt = "n",             # Explicitly remove y-axis
				bty = "n"               # Remove bounding box
			)
			box(bty = "n", col = "white")  # Final clean bounding box removal
			dev.off()
		}
		# reset margins
		par(
			mar = c(5.1, 4.1, 4.1, 2.1)
		)
	}
	## Shapiro test tables for each median metal
	{
		norm_tests$medianmetal <- list()
		for (i in 1:17) {
			x <- metal_points$complete[[13 + i]]
			# if all concentrations identical
			if (length(unique(x)) == 1) {
				norm_tests$medianmetal[[i]] <- data.frame(
					Metal = colnames(metal_points$complete)[13 + i],
					W = NA,
					p.value = NA,
					n = sum(!is.na(x)),
					Note = "Constant values"
				)
				# if concentrations show variance
			} else {
				test <- shapiro.test(x)
				norm_tests$medianmetal[[i]] <- data.frame(
					Metal = colnames(metal_points$complete)[13 + i],
					W = test$statistic,
					p.value = test$p.value,
					n = sum(!is.na(x)),
					Note = ""
				)
			}
		}
		rm(test)
		norm_tests$medianmetal <- do.call(rbind, norm_tests$medianmetal)
		# export table
		write_xlsx(
			x = norm_tests$medianmetal,
			path = "outputs/tables/distest_medianmetal.xlsx"
		)
	}
}
# richness ----
{
	print(shapiro.test(a_div[[5]]$richness))
	hist(a_div[[5]]$richness)
}
# shannon ----
{
	print(shapiro.test(a_div[[5]]$shannon))
	hist(a_div[[5]]$shannon)
}
