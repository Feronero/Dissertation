# Dependencies ----

	{
#		install.packages("vegan")
		library(vegan)
#		install.packages("stringr")
		library(stringr)
#		install.packages("ggplot2")
		library(ggplot2)
#		install.packages("dplyr")
		library(dplyr)
	}

# User-Defined Terms ----
	
	{
		# catchment names in desired order
		catchments <- c(
			"hayle",
			"cober",
			"red",
			"carnon"
		)
		# water sample types of interest
		water_types_of_interest <- c(
			"MINEWATER",
			"GROUNDWATER",
			"RIVER / RUNNING SURFACE WATER",
			"WATER"
		)
		# water quality parameters of interest
		params_of_interest <- c(
			"Copper",
			"Iron",
			"Cadmium"
		)
	}

# Data Import ----

	## Import EA Data ----
	{
		# input oldest year
		oldest <- 2013
		# input newest year
		newest <- 2024
		
		# initialise data storage list
		test_results <- list()
		# step through each year
		for (i in oldest:newest) {
			# sequentially import each dataset
			test_results[[i - (oldest - 1)]] <- read.csv(
				paste0(
					"raw_data/DC-", i, ".csv"
				)
			)
		}
		# rbind each list item together
		test_results <- do.call(
			rbind, test_results
		)
		# create separate year, month, time columns
		test_results[c('year', 'month', 'time')] <- str_split_fixed(
			test_results$sample.sampleDateTime, '-', 3)
		# create year_month column
		test_results$point_year_month <- paste(
			test_results$sample.samplingPoint.notation, test_results$year, test_results$month,
			sep = "_"
		)
		# create parameter_unit column
		test_results$measure_unit <- paste(
			test_results$determinand.definition, test_results$determinand.unit.label,
			sep="_"
		)
		# tidy up measurement units
		test_results$measure_unit <- test_results$measure_unit %>%
			gsub("\\s+", "", .) %>%
			gsub("-", "_", .) %>%
			gsub("/", ".", .)
		# clean up temporary objects
		rm(oldest, newest)
	}

	## Import Sampling Points by Catchment ----
	{
		# initialise storage list
		all_points <- list()
		for (i in 1:length(catchments)) {
			# sequentially import each set of sampling points
			all_points[[i]] <- read.csv(
				paste0(
					"raw_data/sp_", catchments[i], ".csv"
				),
				header = TRUE
			)
			# add a catchment column
			all_points[[i]]$catchment <- rep(
				x = catchments[i],
				length.out = nrow(all_points[[i]])
			)
		}
		# rbind all sampling points together
		all_points <- do.call(
			rbind, all_points
		)
		# remove any sampling points from the wrong area
#		all_points <- subset(
#			x = all_points,
#			subset = all_points$area.label == "Devon and Cornwall"
#		)
	}

# Data Wrangling ----

	## Water Quality Parameters ----
	{
		# initialise parameter storage list
		parameters <- list(
			# extract and store all parameters
			all_params = unique(
				test_results$measure_unit
			),
			selected_params = character()
		)
		# extract and store parameters of interest
		for (i in 1:length(params_of_interest)) {
			parameters$selected_params <- c(
				parameters$selected_params,
				parameters$all_params[
					grep(
						pattern = params_of_interest[i],
						x = parameters$all_params
					)
				]
			)
		}
	}


	## Refine Data ----
	{
		# remove un-needed test_results columns
		filtered_results <- test_results
		filtered_results <- filtered_results[-c(1:9,11,12,14:20)]
		colnames(filtered_results)[2] <- "water_type"
		# keep only river water measurements
		filtered_results <- subset(
			x = filtered_results,
			subset = water_type == "RIVER / RUNNING SURFACE WATER"
		)
		# remove water_type column
		filtered_results <- filtered_results %>% select(-water_type)
		# keep only desired parameter measurements
		filtered_results <- subset(
			x = filtered_results,
			subset = grepl(
				pattern = paste0(
					parameters$selected_params,
					collapse = "|"
				),
				x = filtered_results$measure_unit
			)
		)
		# transpose to wide format
		filtered_results <- reshape(
			data = filtered_results,
			idvar = "point_year_month",
			timevar = "measure_unit",
			direction = "wide"
		)
		# separate point_year_month
		filtered_results <- cbind(
			as.data.frame(
				str_split_fixed(
					string = filtered_results$point_year_month,
					pattern = "_",
					n = 3
				)
			),
			filtered_results
		)
		# name new columns
		colnames(filtered_results)[1] <- "notation"
		colnames(filtered_results)[2] <- "year"
		colnames(filtered_results)[3] <- "month"
		# remove old column
		filtered_results <- filtered_results %>% select(-point_year_month)
	}

	## Extract Sampling Point Metadata ----
	{
		filtered_points <- all_points[,c("notation", "lat", "long", "label", "catchment")]
		filtered_points <- subset(
			x = filtered_points,
			subset = grepl(
				pattern = paste0(
					unique(filtered_results$notation),
					collapse = "|"
				),
				x = filtered_points$notation
			)
		)
	}

	## Merge Filtered Sampling Point Metadata & Test Result Data ----
	{
		filtered_results <- merge(
			x = filtered_results,
			y = filtered_points[, c("notation", "label", "catchment")],
			by = "notation"
		)
	}

# Errors and Questions ----

	# I filtered sampling points by catchment using the URL search in my github readme. Is the method suitable?
	# What parameters did you specify in "some_metals_nutrients" (line 84)
	# error when widening each catchment's table due to duplicate parameter values. Did this occur for you?
		# tested this using your script, appears so
	# site codes in levels(catchment_params_wide$hayle) (line 223 of my script) do not match the string in your script (line 202). Does this indicate a catchment selection error?
	# Some sampling points are from wrong area - remove them after mapping to confirm