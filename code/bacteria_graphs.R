# Alpha-Diversity Plots ----
	
	## box and whisker plots ----
	{
		selection <- "Catchment"
		numeric <- F
		
		# initialise storage list
		merged$alpha.bact_point_copy <- merged$alpha.bact_point
		if (numeric == T) {
			for (i in 1:keyvals$n_levels) {
				# extract factor levels, order them numerically, convert them back to factors
				factor <- factor(
					coda_alpha_copy[[i]][, selection],
					levels = as.factor(
						sort(
							as.numeric(
								levels(
									as.factor(
										coda_alpha_copy[[i]][, selection]
									)
								)
							)
						)
					)
				)
			}
		}
		if (selection == "Month") {
			factor <- droplevels(
				factor(
					x = coda_points[[1]][, selection],
					levels = month.name,
					ordered = TRUE
				)
			)
		} else {
			factor <- coda_points[[1]][, selection]
		}
		# draw box and whisker plots for each taxonomic level
		par(mfrow = c(2, 2))
		for (i in 1:n_levels) {
			for (j in 2:5) {
				boxplot(
					formula = coda_alpha[[i]][, j] ~ factor,
					data = coda_alpha[[i]],
					xlab = selection,
					ylab = colnames(
						coda_alpha_copy[[i]][j]
					)
				)
			}
		}
		# clean up temporary objects
		rm(factor)
	}
	
	## box and whisker plots TEMPORARY FIX ----
	{
		selection <- "Catchment"
		numeric <- F
		
		# initialise storage list
		coda_alpha_copy <- coda_alpha
		if (numeric == T) {
			for (i in 1:n_levels) {
				# extract factor levels, order them numerically, convert them back to factors
				factor <- factor(
					coda_alpha_copy[[i]][, selection],
					levels = as.factor(
						sort(
							as.numeric(
								levels(
									as.factor(
										coda_alpha_copy[[i]][, selection]
									)
								)
							)
						)
					)
				)
			}
		}
		if (selection == "Month") {
			factor <- droplevels(
				factor(
					x = coda_points[[1]][, selection],
					levels = month.name,
					ordered = TRUE
				)
			)
		} else {
			NULL
		}
		# draw box and whisker plots for each taxonomic level
		par(mfrow = c(2, 2))
		for (i in 1:n_levels) {
			for (j in 2:5) {
				boxplot(
					formula = coda_alpha_copy[[i]]$Catchment ~ factor,
					data = coda_alpha[[i]],
					xlab = selection,
					ylab = colnames(
						coda_alpha_copy[[i]][j]
					)
				)
			}
		}
		# clean up temporary objects
		rm(factor)
	}

# Custom nMDS Plots ----

	## plain plot ----
	for (i in 1:n_levels) {
		plot(
			formula = MDS2 ~ MDS1,
			xlim = range(coda_points[[n_levels]][, "MDS1"]),
			ylim = range(coda_points[[n_levels]][, "MDS2"]),
			type = "p",
			data = coda_points[[i]]
		)
	}

	## Colour By Factor ----
	{
		# select grouping factor
		factor <- "Temp"
		# enable user-defined binning - only set true if numeric!
		binning <- T
		# define bin number
		n_bin <- 3
		# gradient or distinct colours?
		gradient_palette <- FALSE
		
		# encode grouping vector based on selection ----
		if (binning == T) {
			# convert data to numeric, divide into bins, convert to factor
			grouping_vector <- as.factor(
				cut(
					x = as.numeric(
						coda_points[[1]][, factor]
					),
					include.lowest = T,
					breaks = n_bin
				)
			)
			# produce readable factor level ranges ----
			range_key <- data.frame(
				value = as.numeric(
					coda_points[[1]][, factor]
				),
				bin = cut(
					x = as.numeric(
						coda_points[[1]][, factor]
					),
					include.lowest = T,
					label = F,
					breaks = n_bin
				),
				notation = grouping_vector
			)
			# order by value
			range_key <- range_key[
				order(
					-range_key$value
				),
			]
			# remove duplicate rows
			range_key <- range_key[
				!duplicated(
					range_key$value
				),
			]
			readable <- character()
			for (i in seq(0, ((2 * n_bin) - 2), 2)) {
				for (j in 1:2) {
					readable[i + j] <- as.character(
						range(
							range_key$value[
								range_key$bin == ((i/2) + 1)
							]
						)[j]
					)	
				}
			}
			legend_text <- paste(
				readable[seq(1, length(readable), by = 2)], 
				readable[seq(2, length(readable), by = 2)], 
				sep = " - "
			)
		} else {
			# order chronologically
			if (factor == "Month") {
				grouping_vector <- droplevels(
					factor(
						x = coda_points[[1]][, factor],
						levels = month.name,
						ordered = TRUE
					)
				)
				legend_text <- levels(grouping_vector)
			} else {
				# convert data directly to factor
				grouping_vector <- as.factor(
					coda_points[[1]][, factor]
				)
				legend_text <- levels(grouping_vector)
			}
		}
		# define grouping colours ----
		if (gradient_palette == TRUE) {
			grouping_colours <- hcl.colors(
				# number of colours scales with factor level
				n = nlevels(grouping_vector),
				palette = "Lajolla"
			)
			# draw plot ----
			par(mfrow = c(1,1))
			for (i in 1:n_levels) {
				plot(
					formula = MDS2 ~ MDS1,
					xlim = range(coda_points[[n_levels]][, "MDS1"]),
					ylim = range(coda_points[[n_levels]][, "MDS2"]),
					col = "black",
					type = "p",
					pch = 21,
					bg = grouping_colours[grouping_vector],
					lwd = 0.1,
					cex = 1.1,
					data = coda_points[[i]]
				)
				# add legend
				legend(
					x = "topleft",
					legend = legend_text,
					col = "black",
					pch = 21,
					pt.bg = grouping_colours,
					pt.lwd = 0.1,
					pt.cex = 1.1,
					title = factor
				)
			}
		} else {
			grouping_colours <- hcl.colors(
				# number of colours scales with factor level
				n = nlevels(grouping_vector),
				palette = "Dark 3"
			)
			# draw plot ----
			par(mfrow = c(1,1))
			for (i in 1:n_levels) {
				plot(
					formula = MDS2 ~ MDS1,
					xlim = range(coda_points[[n_levels]][, "MDS1"]),
					ylim = range(coda_points[[n_levels]][, "MDS2"]),
					col = "black",
					type = "p",
					pch = 21,
					cex = 1.1,
					bg = grouping_colours[grouping_vector],
					lwd = 0.1,
					data = coda_points[[i]]
				)
				# add legend
				legend(
					x = "topleft",
					legend = legend_text,
					col = "black",
					pch = 21,
					pt.cex = 1.1,
					pt.bg = grouping_colours,
					pt.lwd = 0.1,
					title = factor
				)
			}
		}
	}
