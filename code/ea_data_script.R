# Dependencies ----

	{
#		install.packages("vegan")
		library(vegan)
#		install.packages("stringr")
		library(stringr)
#		install.packages("ggplot2")
		library(ggplot2)
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
		# create sepearate year, month, time columns
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

	## Remove Un-Needed Columns ----
	{
		ea_data_trimmed <- ea_data
		# trims identically to tamsyn, except replaces determinand.label with determinand.definition
		ea_data_trimmed <- ea_data_trimmed[
			-c(1,2,4,6,8,9,11,14:17)
		]
	}

	## Import Sampling Points by Catchment ----
	{
		# initialise storage list
		sp_by_catchment <- list()
		for (i in 1:4) {
			# sequentially import each set of sampling points
			sp_by_catchment[[i]] <- read.csv(
				paste0(
					"raw_data/sp_", catchments[i], ".csv"
				),
				header = TRUE
			)
		}
		# name each set by their catchment
		names(sp_by_catchment) <- catchments
	}

# Data Wrangling ----

	## Extract and Select Water Quality Parameters ----
	{
		# initialise parameter storage list
		parameters <- list(
			# extract and store all parameters
			all_params = unique(
				ea_data_trimmed$determinand.definition
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
		# clean up temporary objects
		rm(params_of_interest)
	}

	## Subset EA Data by Water Type ----
	{
		# initialise storage list
		water_types <- list()
		# extract and store samples from each water type of interest
		for (i in 1:length(water_types_of_interest)) {
			water_types[[i]] <- ea_data_trimmed[
				grep(
					pattern = water_types_of_interest[i],
					x = ea_data_trimmed$sample.sampledMaterialType.label
				),
			]
			# remove redundant "water type" column
			water_types[[i]] <- water_types[[i]][, -6]
		}
		# name each water type
		names(water_types) <- water_types_of_interest
	}

	## Subset River Water by Catchment ----
	{
		# initialise storage list
		river_by_catchment <- list()
		for (i in 1:length(catchments)) {
			river_by_catchment[[i]] <- subset(
				x = water_types$`RIVER / RUNNING SURFACE WATER`,
				subset = grepl(
					pattern = paste(
						sp_by_catchment[[i]]$notation,
						collapse = "|"
					),
					x = water_types$`RIVER / RUNNING SURFACE WATER`$sample.samplingPoint.notation
				)
			)
		}
		names(river_by_catchment) <- catchments
	}
	
	## Subset Catchment by selected_parameters ----
	{
		catchment_and_params <- list()
		for (i in 1:length(catchments)) {
			catchment_and_params[[i]] <- subset(
				x = river_by_catchment[[i]],
				subset = grepl(
					pattern = paste(
						parameters$selected_params,
						collapse = "|"
					),
					x = river_by_catchment[[i]]$determinand.definition
				)
			)
		}
		# appropriately name list items
		names(catchment_and_params) <- catchments
		# remove un-need columns
		for(i in 1:length(catchments)) {
			catchment_and_params[[i]] <- subset(
				x = catchment_and_params[[i]],
				select = c(
					result,
					Sample_Type_Time,
					Measure_Unit
				)
			)
		}
	}
	
	## Transpose Catchment Data to Wide Format ----
	{
		catchment_params_wide <- list()
		for (i in 1:length(catchments)) {
			catchment_params_wide[[i]] <- reshape(
				data = catchment_and_params[[i]],
				idvar = "Sample_Type_Time",
				timevar = "Measure_Unit",
				direction = "wide"
			)
		}
		names(catchment_params_wide) <- catchments
	}

	
new_cols <- as.data.frame(str_split_fixed(catchment_params_wide$hayle$Sample_Type_Time, "_", 2))
wide_hayle$Sample <-paste(new_cols$V1, new_cols$V2, sep="_")
wide_hayle$Site <- new_cols$V1
wide_hayle$Time <- new_cols$V2

new_cols1 <- as.data.frame(str_split_fixed(wide_red$Sample_Type_Time, "_", 2))
wide_red$Sample <-paste(new_cols1$V1, new_cols1$V2, sep="_")
wide_red$Site <- new_cols1$V1
wide_red$Time <- new_cols1$V2

new_cols2 <- as.data.frame(str_split_fixed(wide_carnon$Sample_Type_Time, "_", 2))
wide_carnon$Sample <-paste(new_cols2$V1, new_cols2$V2, sep="_")
wide_carnon$Site <- new_cols2$V1
wide_carnon$Time <- new_cols2$V2

new_cols3 <- as.data.frame(str_split_fixed(wide_cober$Sample_Type_Time, "_", 2))
wide_cober$Sample <-paste(new_cols3$V1, new_cols3$V2, sep="_")
wide_cober$Site <- new_cols3$V1
wide_cober$Time <- new_cols3$V2

# Calculate Linear Regressions ----

	{
		
	}

# Errors and Questions ----

	# I filtered sampling points by catchment using the URL search in my github readme. Is the method suitable?
	# What parameters did you specify in "some_metals_nutrients" (line 84)
	# error when widening each catchment's table due to duplicate parameter values. Did this occur for you?