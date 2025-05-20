# Merges EA metal data with metadata from the appropriate sampling point to link metal amounts with location, altitude, distance to mouth etc. See "merged$ea_data.points"

# 0.0 Dependencies ----

{
#	install.packages("vegan")
	library(vegan)
#	install.packages("stringr")
	library(stringr)
#	install.packages("ggplot2")
	library(ggplot2)
#	install.packages("dplyr")
	library(dplyr)
}

# 2.1 Key Terms ----
	
{
	keyvals <- list(
		# catchment names in desired order
		catchments = c(
			"hayle",
			"cober",
			"red",
			"carnon"
		),
		# water sample types of interest
		water_types = c(
			"MINEWATER",
			"GROUNDWATER",
			"RIVER / RUNNING SURFACE WATER",
			"WATER"
		),
		# water quality parameters of interest
		parameters = c(
			"Copper",
			"Iron",
			"Cadmium",
			"Arsenic"
		),
		years = c(
			2013,
			2024
		)
	)
}

# Data Import ----

	## Import EA Data ----
		## creates ea data list, and inserts raw EA data into $all
		## adds a few columns to aid later transposition to wide format
	{
		# initialise data storage list
		ea_data <- list()
		# step through each year
		for (i in keyvals$years[1]:keyvals$years[2]) {
			# sequentially import each dataset
			ea_data$all[[i - (keyvals$years[1] - 1)]] <- read.csv(
				paste0(
					"raw_data/EA_water_data/DC-", i, ".csv"
				)
			)
		}
		# rbind each list item together
		ea_data$all <- do.call(
			rbind, ea_data$all
		)
		# create separate year, month, time columns
		ea_data$all[c('year', 'month', 'time')] <- str_split_fixed(
			ea_data$all$sample.sampleDateTime, '-', 3)
		# create "year_month" column
		ea_data$all$point_year_month <- paste(
			ea_data$all$sample.samplingPoint.notation, ea_data$all$year, ea_data$all$month,
			sep = "_"
		)
		# create "measure_unit" column
		ea_data$all$measure_unit <- paste(
			ea_data$all$determinand.definition, ea_data$all$determinand.unit.label,
			sep="_"
		)
		# tidy up "measurement units" contents
		ea_data$all$measure_unit <- ea_data$all$measure_unit %>%
			gsub("\\s+", "", .) %>%
			gsub("-", "_", .) %>%
			gsub("/", ".", .)
	}

	## Import Sampling Points by Catchment ----
		## creates points list, inserts all relevant points into $all
		## labels each point with its catchment name
	{
		# initialise storage list
		ea_points <- list()
		# fill "all"
		for (i in 1:length(keyvals$catchments)) {
			# sequentially import each set of sampling points
			ea_points$all[[i]] <- read.csv(
				paste0(
					"raw_data/sampling_points/sp_", keyvals$catchments[i], ".csv"
				),
				header = TRUE
			)
			# add a catchment column
			ea_points$all[[i]]$catchment <- rep(
				x = keyvals$catchments[i],
				length.out = nrow(ea_points$all[[i]])
			)
		}
		# rbind all sampling points together
		ea_points$all <- do.call(
			rbind, ea_points$all
		)
		# set "catchment" column to factor
		ea_points$all$catchment <- as.factor(
			ea_points$all$catchment
		)
	}

# Data Wrangling ----

	## Filter EA Water Quality Parameters ----
	{
		# initialise parameter storage list
		measure_units <- list()
		# extract and store all EA water quality parameters
		measure_units$all = unique(
			ea_data$all$measure_unit
		)
		# extract and store parameters of interest
		for (i in 1:length(keyvals$parameters)) {
			measure_units$filtered <- c(
				measure_units$filtered,
				measure_units$all[
					grep(
						pattern = keyvals$parameters[i],
						x = measure_units$all
					)
				]
			)
		}
	}


	## Filter EA Data by Selected Parameters ----
	{
		# prepare a $filtered section
		ea_data$filtered <- ea_data$all
		# remove un-needed ea_data$filtered columns left over from ea_data$all
		ea_data$filtered <- ea_data$filtered[, c(21,10,22)]
#		colnames(ea_data$filtered)[2] <- "water_type"
#		# keep only river water measurements
#		ea_data$filtered <- subset(
#			x = ea_data$filtered,
#			subset = water_type == "RIVER / RUNNING SURFACE WATER"
#		)
#		# remove water_type column
#		ea_data$filtered <- ea_data$filtered %>% select(-water_type)
		# keep only desired parameter measurements
		ea_data$filtered <- subset(
			x = ea_data$filtered,
			subset = grepl(
				pattern = paste0(
					measure_units$filtered,
					collapse = "|"
				),
				x = ea_data$filtered$measure_unit
			)
		)
	}
	## Transpose EA Data to Wide Format ----
	{
		# transpose to wide format
		ea_data$wide <- reshape(
			data = ea_data$filtered,
			idvar = "point_year_month",
			timevar = "measure_unit",
			direction = "wide"
		)
		# separate point_year_month
		ea_data$wide <- cbind(
			as.data.frame(
				str_split_fixed(
					string = ea_data$wide$point_year_month,
					pattern = "_",
					n = 3
				)
			),
			ea_data$wide
		)
		# name new columns
		colnames(ea_data$wide)[1] <- "notation"
		colnames(ea_data$wide)[2] <- "year"
		colnames(ea_data$wide)[3] <- "month"
		# remove now-redundant point_year_month column
		ea_data$wide <- ea_data$wide[, -4]
	}

	## Filter Sampling Point Metadata ----
	{
		# creates $filtered item in points list
		ea_points$filtered <- ea_points$all
		# keeps all points whose notation code is shared with the filtered ea_data
		ea_points$filtered <- subset(
			x = ea_points$filtered,
			subset = grepl(
				pattern = paste0(
					unique(ea_data$wide$notation),
					collapse = "|"
				),
				x = ea_points$filtered$notation
			)
		)
		# coerces the catchment column to type factor
		ea_points$filtered$catchment <- as.factor(
			ea_points$filtered$catchment
		)
	}

	## Merge Filtered Sampling Point Metadata & Test Result Data ----
	{
		merged <- list()
		merged$ea_data.points <- merge(
			x = ea_data$wide,
			y = ea_points$filtered,
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



# Operate under assumption that ea data and bacterial data are paired? YES
# are the pollutants essential metals? Essential metals may be more easily bio-absorbed due to existing metabolic pathways
# look for clustering of certain metals and env. parameters
# if time, investigate seasonal effects?


# add the following data to ea_data$filtered:
	# altitude
	# distance to mouth
	# order from source to mouth