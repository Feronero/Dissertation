{
	## Metadata for all bacteria samples
		bact_data <- list()
		bact_data$all <- read.csv(
			"raw_data/sample_metadata_blanks_removed.csv"
		)
		names(bact_data$all)[15] <- "DTM.km"
	## Metadata for all bacteria points
		bact_points <- list()
		bact_points$all <- bact_data$all[
			!duplicated(bact_data$all$Site),
			c(3, 4, 5, 12, 13, 14, 15, 16, 20)
		]
	## PCR Data
		PCR_data <- list()
		for (i in 1:keyvals$n_levels) {
			PCR_data[[i]] <- read.csv(
				paste0(
					"raw_data/PCR/", keyvals$PCR_levels[[i]], ".csv"
				)
			)
			colnames(PCR_data[[i]])[1] <- "SampleID"
		}
	## qPCR Data
		qPCR_data <- list()
		for (i in 1:7) {
			qPCR_data$all[[i]] <- read.csv(
				paste0(
					"raw_data/qPCR_refined/plate_", i, ".csv"
				)
			)[, 1:2]
		qPCR_data$all <- do.call(
			rbind, qPCR_data$all
		)
		# distill out the mean of each sample pair
			qPCR_data$meaned <- qPCR_data$all %>%
				group_by(Sample) %>%
				summarise(
					Mean_Cq = mean(Cq, na.rm = TRUE),
					SD_Cq = sd(Cq, na.rm = TRUE)
				)
	}
	
	## Alpha Diversity
	{
		a_div <- list()
		for (i in 1:keyvals$n_levels) {
			a_div[[i]] <- data.frame(
				SampleID = PCR_data[[i]][,"SampleID"]
			)
		}
		# species richness
		for (i in 1:keyvals$n_levels) {
			a_div[[i]][,"richness"] <- specnumber(
				PCR_data[[i]][, -1]
			)
		}
		# all other a_indices
		for (j in keyvals$a_indices[-1]) {
			for (i in 1:keyvals$n_levels) {
				a_div[[i]][, j] <- diversity(
					PCR_data[[i]][, -1],
					index = j
				)
			}
		}
	}
}
