{ # Import Raw Bacteria Data ----

	# Import metadata for all bacteria samples
		bact_meta <- list()
		bact_meta <- read.csv(
			"raw_data/sample_metadata_blanks_removed.csv"
		)
		names(bact_meta)[15] <- "DTM.km"
		bact_meta$Year <- as.numeric(format(as.Date(bact_meta$Date, format = "%d/%m/%Y"), "%Y"))
		bact_meta$Month <- format(as.Date(bact_meta$Date, format = "%d/%m/%Y"), "%m")
		bact_meta <- bact_meta %>%
			select(
				-Grid_ref,
				-X,
				-Y,
				-Latitude,
				-Longitude,
				-Time,
				-Altitude,
				-EC,
				-Site_code,
				-Site_code_2,
				-Site_time_code,
				-Poll_class
			) %>%
			relocate(
				Year,
				.after = "Location"
			)
	# Import PCR Data
		PCR_data <- list()
		for (i in 1:keyvals$n_levels) {
			PCR_data[[i]] <- read.csv(
				paste0(
					"raw_data/PCR/", keyvals$PCR_levels[[i]], ".csv"
				)
			)
			colnames(PCR_data[[i]])[1] <- "SampleID"
		}
	# Import qPCR Data
		qPCR_data <- list()
		for (i in 1:7) {
			qPCR_data[[i]] <- read.csv(
				paste0(
					"raw_data/qPCR_refined/plate_", i, ".csv"
				)
			)[, 1:2]
		}
		qPCR_data <- do.call(
			rbind, qPCR_data
		)
}

{ # Calculate Mean qPCR Cq Values (from the two technical replicates each) ----

	qPCR_data <- qPCR_data %>%
		group_by(Sample) %>%
		summarise(
			Mean_Cq = mean(Cq, na.rm = TRUE),
			SD_Cq = sd(Cq, na.rm = TRUE)
		)
}

{ # Calculate Alpha Diversity ----

	# create alpha_diversity storage list
		a_div <- list()
		for (i in 1:keyvals$n_levels) {
			a_div[[i]] <- data.frame(
				SampleID = PCR_data[[i]][,"SampleID"]
			)
		}
	# calculate species richness
		for (i in 1:keyvals$n_levels) {
			a_div[[i]][,"richness"] <- specnumber(
				PCR_data[[i]][, -1]
			)
		}
	# calculate all other a_indices
		for (j in keyvals$a_indices[-1]) {
			for (i in 1:keyvals$n_levels) {
				a_div[[i]][, j] <- diversity(
					PCR_data[[i]][, -1],
					index = j
				)
			}
		}
}

{ # Merge all data from each bacteria sample into single dataframe ----
	
	# create list
		bact_data <- list()
	# merge bact_meta and qPCR data for each tax. level
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- merge(
				x = bact_meta,
				y = qPCR_data,
				by = "Sample",
				all.x = TRUE
			)
		}
	# add alpha_div data
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- merge(
				x = bact_data[[i]],
				y = a_div[[i]],
				by = "SampleID"
			)
		}
	# add PCR data
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- merge(
				x = bact_data[[i]],
				y = PCR_data[[i]],
				by = "SampleID"
			)
		}
	# add "season" factor column
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- bact_data[[i]] %>%
				mutate(
					Season = case_when(
						Month %in% "04" ~ "Spring",
						Month %in% "07" ~ "Summer",
						Month %in% "10" ~ "Autumn",
						Month %in% "01" ~ "Winter"
					),
					Season = factor(
						Season,
						levels = c(
							"Spring", "Summer", "Autumn", "Winter"
						)
					)
				) %>%
				relocate(
					Season,
					.after = 5
				)
		}
	# export the level 4 bact_data as a .csv
		write_csv(
			x = bact_data[[1]],
			file = "outputs/tables/bact_data.csv"
		)
	# remove individual dataframes
		rm(
			a_div, bact_meta, PCR_data, qPCR_data
		)
	
}

{ # Check subsampling depth
	f = 1
	rowSums(bact_data[[f]][1:nrow(bact_data[[f]]), 21:ncol(bact_data[[f]])])
}
