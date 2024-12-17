# Dependencies ----

	{
#		install.packages("vegan")
		library(vegan)
#		install.packages("stringr")
		library(stringr)
#		install.packages("ggplot2")
		library(ggplot2)
	}

# Data Import ----

	# Import Sample Data ----
	{
		# import sample metadata
		sample_metadata <- read.csv(
			file = "raw_data/sample_metadata_blanks_removed.csv",
			stringsAsFactors = TRUE
		)
		# import ASV data across all taxonomic levels
		asv_table <- read.csv(
			"raw_data/asv_table.csv"
		)
		# import ASV data filtered by each taxonomic level, and bundle them into a list
		ASVs <- list(
			level_2 = read.csv("raw_data/level_2.csv"),
			level_3 = read.csv("raw_data/level_3.csv"),
			level_4 = read.csv("raw_data/level_4.csv"),
			level_5 = read.csv("raw_data/level_5.csv"),
			level_6 = read.csv("raw_data/level_6.csv")
		)
		# rename the index column to match "sample_metadata"
		for (i in 1:length(ASVs)) {
			colnames(ASVs[[i]])[1] <- "SampleID"
		}
		# check output
		print(ASVs[[1]][1:4,1:4])
	}

	# Import EA Data ----
	{
		# input oldest year
		oldest <- 2013
		# input newest year
		newest <- 2024
		
		# initialise data storage list
		ea_data <- list()
		# step through each year
		for (i in oldest:newest) {
			# sequentially name and import each dataset
			ea_data[[i]] <- read.csv(
				paste0(
					"raw_data/DC-", i, ".csv"
				)
			)
		}
		# use do.call to rbind each list item together
		ea_data <- do.call(
			rbind, ea_data
		)
		# create seperate year, month, time columns
		ea_data[c('Year', 'Month', 'Time')] <- str_split_fixed(
			ea_data$sample.sampleDateTime, '-', 3)
		# create year_month column
		ea_data$Sample_Type_Time <- paste(
			ea_data$sample.samplingPoint.notation, ea_data$Year, ea_data$Month,
			sep = "_"
		)
		# create parameter_unit column
		ea_data$Measure_Unit <- paste(
			ea_data$determinand.label, ea_data$determinand.unit.label,
			sep="_"
		)
		# clean up temporary objects
		rm(oldest, newest)
	}

	# Import Sampling Points by Catchment ----
	{
		sp_by_catchment <- list(
			sp_hayle = read.csv(
				"raw_data/sp_hayle.csv",
				header = TRUE
			),
			sp_cober = read.csv(
				"raw_data/sp_cober.csv",
				header = TRUE
			),
			sp_red = read.csv(
				"raw_data/sp_red.csv",
				header = TRUE
			),
			sp_carnon = read.csv(
				"raw_data/sp_carnon.csv",
				header = TRUE
			)	
		)
	}

# Data Wrangling ----

	## Order Metadata Factor Levels ----
	{
		# order "Month" column
		sample_metadata$Month <- droplevels(
			factor(
				x = as.factor(
					sample_metadata$Month
				),
				levels = month.name,
				ordered = TRUE
			)
		)
	}

	## Extract and Select Water Quality Parameters ----
	{
		# extract all parameters
		all_params = data.frame(
			all_params = unique(
				ea_data$determinand.definition
			)
		)
		# define selected_parameters
		selected_params = data.frame(
			selected_params = c(
				# copper
				all_params$all_params[
					grep(
						pattern = "Copper",
						x = all_params$all_params
					)
				],
				# iron
				all_params$all_params[
					grep(
						pattern = "Iron",
						x = all_params$all_params
					)
				],
				# cadmium
				all_params$all_params[
					grep(
						pattern = "Cadmium",
						x = all_params$all_params
					)
				]
			)
		)
	}

	## Subset EA Data by Water Type ----
	{
		# subset by row where sample type == MINEWATER
		MINEWATER <- ea_data[
			grep(
				pattern = "MINEWATER",
				x = ea_data$sample.sampledMaterialType.label
			),
		]
		# subset by row where sample type == GROUNDWATER
		GROUNDWATER <- ea_data[
			grep(
				pattern = "GROUNDWATER",
				x = ea_data$sample.sampledMaterialType.label
			),
		]
		# subset by row where sample type == SURFACE / RUNNING SURFACE WATER
		RIVER <- ea_data[
			grep(
				pattern = "RIVER / RUNNING SURFACE WATER",
				x = ea_data$sample.sampledMaterialType.label
			),
		]
		# subset by row where sample type == WATER
		ALL_WATER <- ea_data[
			grep(
				pattern = "WATER",
				x = ea_data$sample.sampledMaterialType.label
			),
		]
	}

	## Subset Water Type by Catchment ----
	
		### subset RIVER water by catchment ----
		
			# fetches sampling point IDs from a specific catchment
			# joins fetched IDs into a single vector
			# marks each join with "|", which acts as an OR operator
			# searches RIVER sampling points for rows with matching IDs
			# saves matching rows to a catchment-specific object
			# RIVER is now filtered by catchment

		{
			# search RIVER for sampling point IDs in the "sp_hayle" data frame
			RIVER_hayle <- RIVER[
				grepl(
					pattern = paste(
						sp_by_catchment[[1]]$notation,
						collapse = "|"
					),
					x = RIVER$sample.samplingPoint.notation
				),
			]
			# search RIVER for sampling point IDs in the "sp_cober" data frame
			RIVER_cober <- RIVER[
				grepl(
					pattern = paste(
						sp_by_catchment[[2]]$northing,
						collapse = "|"
					),
					x = RIVER$sample.samplingPoint.notation
				),
			]
			# search RIVER for sampling point IDs in the "sp_red" data frame
			RIVER_red <- RIVER[
				grepl(
					pattern = paste(
						sp_by_catchment[[3]]$notation,
						collapse = "|"
					),
					x = RIVER$sample.samplingPoint.notation
				),
			]
			# search RIVER for sampling point IDs in the "sp_carnon" data frame
			RIVER_carnon <- RIVER[
				grepl(
					pattern = paste(
						sp_by_catchment[[4]]$notation,
						collapse = "|"
					),
					x = RIVER$sample.samplingPoint.notation
				),
			]
		}
	
		# subset catchment by selected_parameters
		{
			subset_hayle <- RIVER_hayle[
				grepl(
					pattern = paste(
						selected_params,
						collapse = "|"
					),
					x = RIVER_hayle$determinand.definition
				),
			]
			subset_cober <- RIVER_cober[
				grepl(
					pattern = paste(
						selected_params,
						collapse = "|"
					),
					x = RIVER_cober$determinand.definition
				),
			]
			subset_red <- RIVER_red[
				grepl(
					pattern = paste(
						selected_params,
						collapse = "|"
					),
					x = RIVER_red$determinand.definition
				),
			]
			subset_carnon <- RIVER_carnon[
				grepl(
					pattern = paste(
						selected_params,
						collapse = "|"
					),
					x = RIVER_carnon$determinand.definition
				),
			]
		}
	
	## Transpose Catchment Data to Wide Format ----
	{
		subset_hayle_wide <- reshape(
			subset_hayle,
			idvar = "Sample_Type_Time",
			timevar = "Measure_Unit",
			direction = "wide"
		)
		subset_cober_wide <- reshape(
			subset_cober,
			idvar = "Sample_Type_Time",
			timevar = "Measure_Unit",
			direction = "wide"
		)
		subset_red_wide <- reshape(
			subset_red,
			idvar = "Sample_Type_Time",
			timevar = "Measure_Unit",
			direction = "wide"
		)
		subset_carnon_wide <- reshape(
			subset_carnon,
			idvar = "Sample_Type_Time",
			timevar = "Measure_Unit",
			direction = "wide"
		)
	}

# Calculate Linear Regressions ----

	{
		
	}

# Errors and Questions ----

	# I filtered sampling points by catchment using the URL search in my github readme. Is the method suitable?
	# What parameters did you specify in "some_metals_nutrients" (line 84)
	# error when widening each catchment's table due to duplicate parameter values. Did this occur for you?