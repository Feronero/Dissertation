# Metal Sampling Points ----
{
	metal_points <- list()
	# all points
	for (i in 1:length(keyvals$catchments)) {
		metal_points$all[[i]] <- read.csv(
			paste0(
				"raw_data/sampling_points/sp_", keyvals$catchments[i], ".csv"
			),
			header = TRUE
		)
	}
	metal_points$all <- do.call(
		rbind, metal_points$all
	)
	metal_points$all$catchment <- as.factor(
		metal_points$all$catchment
	)
	# metal points with a bacterial point twin
	metal_points$paired <- metal_points$all[
		!is.na(metal_points$all$eDNA_point)
		,
	]
}

# Metal Data ----
{
	## 1. Import data
	{
		metal_data <- list()
		for (i in keyvals$years[1]:keyvals$years[2]) {
			metal_data$all[[i - (keyvals$years[1] - 1)]] <- read.csv(
				paste0(
					"raw_data/EA_water_data/DC-", i, ".csv"
				)
			)
		}
		# rbind each list item together
		metal_data$all <- do.call(
			rbind, metal_data$all
		)
	}
	## 2. Remove un-paired points
	{	
		# remove results from un-paired metal points
		metal_data$paired <- subset(
			x = metal_data$all,
			subset = sample.samplingPoint.notation %in% metal_points$paired$notation
		)
	}
	## 3. Tidy up columns
	{
		# create separate year, month, time columns
		metal_data$paired[c('year', 'month', 'time')] <- str_split_fixed(
			metal_data$paired$sample.sampleDateTime, '-', 3)
		# create "year_month" column
		metal_data$paired$point_year_month <- paste(
			metal_data$paired$sample.samplingPoint.notation, metal_data$paired$year, metal_data$paired$month,
			sep = "_"
		)
		# create "measure_unit" column
		metal_data$paired$measure_unit <- paste(
			metal_data$paired$determinand.definition, metal_data$paired$determinand.unit.label,
			sep="_"
		)
		# tidy up "measurement units" contents
		metal_data$paired$measure_unit <- metal_data$paired$measure_unit %>%
			gsub("\\s+", "", .) %>%
			gsub("-", "_", .) %>%
			gsub("/", ".", .)
	}
	## 4. Filter out non-metal results
	{	
		# remove unnecessary results
		metal_data$filtered <- subset(
			x = metal_data$paired,
			subset = determinand.definition %in% keyvals$parameters
		)
		# remove unnecessary columns
		metal_data$filtered <- metal_data$filtered[, c(21,10,22)]
	}
	## 5. Create list of each metal and its measure unit for reference
	{
		metal_data$reference <- data.frame(
			metals_and_units = unique(
				metal_data$filtered$measure_unit
			)
		)
	}
	## 6. Transpose to wide format
	{	
		metal_data$wide <- reshape(
			data = metal_data$filtered,
			idvar = "point_year_month",
			timevar = "measure_unit",
			direction = "wide"
		)
		# separate point_year_month
		metal_data$wide <- cbind(
			as.data.frame(
				str_split_fixed(
					string = metal_data$wide$point_year_month,
					pattern = "_",
					n = 3
				)
			),
			metal_data$wide
		)
	}
	## 7. Tidy up wide column names
	{
		# name new columns
		colnames(metal_data$wide)[1] <- "notation"
		colnames(metal_data$wide)[2] <- "year"
		colnames(metal_data$wide)[3] <- "month"
		# remove now-redundant point_year_month column
		metal_data$wide <- metal_data$wide[, -4]
		# Loop through desired names and find matching columns
		for (i in keyvals$metals) {
			# Find column index where the name appears (partial match)
			match_index <- grep(
				i,
				colnames(metal_data$wide),
				ignore.case = TRUE
			)
			# If a match is found, rename the column
			if (length(match_index) == 1) {
				colnames(metal_data$wide)[match_index] <- i
			}
		}
	}
	## 8. Convert all metals to μg/L
	{
		for (i in keyvals$milligrams) {
			metal_data$wide[[i]] <- metal_data$wide[[i]] * 1000
		}
	}
	## 9. Add sampling point information
	{
		metal_data$wide <- merge(
			x = metal_data$wide,
			y = metal_points$paired,
			by = "notation"
		)
	}
	## 11. Calculate median median levels per metal point
		## normalise for biological impact if able
	{
		# Create list to store medians and interquartile ranges
		store <- list()
		store$median <- as.data.frame(
			matrix(
				data = NA,
				nrow = nrow(metal_points$paired),
				ncol = length(keyvals$metals)
			)
		)
		store$iqr <- as.data.frame(
			matrix(
				data = NA,
				nrow = nrow(metal_points$paired),
				ncol = length(keyvals$metals)
			)
		)
		names(store$median) <- keyvals$metals
		names(store$iqr) <- keyvals$metals
		
		for (i in seq_along(metal_points$paired$notation)) {
			for (j in keyvals$metals) {
				# Subset data for current sampling point and metal
				subset_data <- subset(
					x = metal_data$wide,
					subset = notation == metal_points$paired$notation[i]
				)[[j]]
				
				# If no data, output NA. Otherwise, compute median and IQR
				if (length(subset_data) == 0 || all(is.na(subset_data))) {
					store$median[[j]][i] <- NA
					store$iqr[[j]][i] <- NA
				} else {
					store$median[[j]][i] <- median(subset_data, na.rm = TRUE)
					store$iqr[[j]][i] <- IQR(subset_data, na.rm = TRUE)
				}
			}
		}
		
		# Assign final column names and append to metal_points
		names(store$median) <- paste0(keyvals$metals, "_median.μgL")
		names(store$iqr) <- paste0(keyvals$metals, "_iqr.μgL")
		metal_points$paired <- cbind(metal_points$paired, store$median, store$iqr)
		# clean up temporary objects
		rm(store, subset_data)
	}
}

# Archive ----
# 
# 	## 12. Classify total metal amounts
# 	{
# 		# 1. Select complete "_median" columns
# 		metal_cols <- grep("_median.μgL$", names(metal_points$paired), value = TRUE)
# 		complete_metal_cols <- metal_cols[colSums(is.na(metal_points$paired[, metal_cols])) == 0]
# 		
# 		# 2. Calculate total median for each sampling point in mg/L
# 		metal_points$paired$Total_Median.mgL <- (rowSums(metal_points$paired[, complete_metal_cols])) / 1000
# 		
# 		
# 		# 3. Classify into quartiles
# 		quartiles <- quantile(
# 			metal_points$paired$Total_Median,
# 			probs = seq(0, 1, 0.25),
# 			na.rm = TRUE
# 		)
# 		metal_points$paired$Quartile_Class <- cut(
# 			metal_points$paired$Total_Median,
# 			breaks = quartiles,
# 			labels = c("Q1", "Q2", "Q3", "Q4"),
# 			include.lowest = TRUE
# 		)
# 	}
# 	## 12. Classify individual metal amounts
# 	{
# 		# Identify the metal mean columns (ending with "_mean.μgl")
# 		metal_mean_cols <- grep(
# 			"_mean\\.μgL$", names(metal_points$paired), value = TRUE
# 		)
# 		
# 		# Create quartile columns for each metal
# 		for (i in metal_mean_cols) {
# 			# Create new column name
# 			quartile_col <- gsub("_mean\\.μgl$", "_quartile", i)
# 			
# 			# Calculate quartile breaks
# 			concentrations <- metal_points$paired[[i]]
# 			quartile_breaks <- quantile(
# 				concentrations, 
# 				probs = c(0, 0.25, 0.5, 0.75, 1), 
# 				na.rm = TRUE
# 			)
# 			
# 			# If all values identical, classify all as Q4
# 			if (length(unique(quartile_breaks)) == 1) {
# 				metal_points$paired[[quartile_col]] <- "Q4"
# 			} else {
# 				# Categorize into quartiles
# 				metal_points$paired[[quartile_col]] <- cut(
# 					concentrations,
# 					breaks = quartile_breaks,
# 					include.lowest = TRUE,
# 					labels = c("Q1", "Q2", "Q3", "Q4")
# 				)
# 			}
# 		}
# 		
# 		# Convert quartile columns to factors (optional)
# 		quartile_cols <- grep(
# 			"_quartile$", names(metal_points$paired), value = TRUE
# 		)
# 		metal_points$paired[quartile_cols] <- lapply(
# 			metal_points$paired[quartile_cols], factor
# 		)
# 	}
