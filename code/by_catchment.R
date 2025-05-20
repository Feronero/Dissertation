# master list of all bacterial PCR samples with their sampling point data
{
	merged$PCR.points <- list()
	for (i in 1:keyvals$n_levels) {
		merged$PCR.points[[i]] <- merge(
			x = bact_data$ASVs[[i]],
			y = bact_points$all,
			by = "SampleID"
		)
	}
}

# Make Storage lists for alpha and beta div by catchment ----
{
	# Hayle
		## all hayle points
		hayle <- list()
		for (i in 1:keyvals$n_levels) {
			hayle$samples[[i]] <- subset(
				x = merged$PCR.points[[i]],
				subset = merged$PCR.points[[i]]$Catchment == "Hayle"
			)
		}
		## number of samples from hayle
		hayle$count <- nrow(hayle$samples[[1]])
		## prepare storage list for hayle alpha diversity
		hayle$a_indices <- list()
		for (i in 1:keyvals$n_levels) {
			hayle$a_indices[[i]] <- data.frame(
				SampleID = hayle$samples[[1]]$SampleID,
				Richness = rep(0, hayle$count),
				Shannon = rep(0, hayle$count),
				Simpson = rep(0, hayle$count),
				Invsimpson = rep(0, hayle$count)
			)
		}
	# Cober
		## all cober points
		cober <- list()
		for (i in 1:keyvals$n_levels) {
			cober$samples[[i]] <- subset(
				x = merged$PCR.points[[i]],
				subset = merged$PCR.points[[i]]$Catchment == "Hayle"
			)
		}
		## number of samples from cober
		cober$count <- nrow(cober$samples[[1]])
		## prepare storage list for cober alpha diversity
		cober$a_indices <- list()
		for (i in 1:keyvals$n_levels) {
			cober$a_indices[[i]] <- data.frame(
				SampleID = cober$samples[[1]]$SampleID,
				Richness = rep(0, cober$count),
				Shannon = rep(0, cober$count),
				Simpson = rep(0, cober$count),
				Invsimpson = rep(0, cober$count)
			)
		}
	# Carnon
		## all carnon points
		carnon <- list()
		for (i in 1:keyvals$n_levels) {
			carnon$samples[[i]] <- subset(
				x = merged$PCR.points[[i]],
				subset = merged$PCR.points[[i]]$Catchment == "Hayle"
			)
		}
		## number of samples from carnon
		carnon$count <- nrow(carnon$samples[[1]])
		## prepare storage list for carnon alpha diversity
		carnon$a_indices <- list()
		for (i in 1:keyvals$n_levels) {
			carnon$a_indices[[i]] <- data.frame(
				SampleID = carnon$samples[[1]]$SampleID,
				Richness = rep(0, carnon$count),
				Shannon = rep(0, carnon$count),
				Simpson = rep(0, carnon$count),
				Invsimpson = rep(0, carnon$count)
			)
		}

	# red
		## all red points
		red <- list()
		for (i in 1:keyvals$n_levels) {
			red$samples[[i]] <- subset(
				x = merged$PCR.points[[i]],
				subset = merged$PCR.points[[i]]$Catchment == "Hayle"
			)
		}
		## number of samples from red
		red$count <- nrow(red$samples[[1]])
		## prepare storage list for red alpha diversity
		red$a_indices <- list()
		for (i in 1:keyvals$n_levels) {
			red$a_indices[[i]] <- data.frame(
				SampleID = red$samples[[1]]$SampleID,
				Richness = rep(0, red$count),
				Shannon = rep(0, red$count),
				Simpson = rep(0, red$count),
				Invsimpson = rep(0, red$count)
			)
		}
}

# Calculate Per-Catchment Alpha-Diversity Indices ----
{
	## hayle
	
		# species richness
		for (i in 1:keyvals$n_levels) {
			hayle$a_indices[[i]][, 2] <- specnumber(
				hayle$samples[[i]][, 2:(ncol(hayle$samples[[i]]) - 22)]
			)
		}
		# Shannon diversity
		for (i in 1:keyvals$n_levels) {
			hayle$a_indices[[i]][, 3] <- diversity(
				hayle$samples[[i]][, 2:(ncol(hayle$samples[[i]]) - 22)],
				index = "shannon"
			)
		}
		# Simpson diversity
		for (i in 1:keyvals$n_levels) {
			hayle$a_indices[[i]][, 4] <- diversity(
				hayle$samples[[i]][, 2:(ncol(hayle$samples[[i]]) - 22)],
				index = "simpson"
			)
		}
		# inverse Simpson
		for (i in 1:keyvals$n_levels) {
			hayle$a_indices[[i]][, 5] <- diversity(
				hayle$samples[[i]][, 2:(ncol(hayle$samples[[i]]) - 22)],
				index = "invsimpson"
			)
		}
	
	## cober
	
		# species richness
		for (i in 1:keyvals$n_levels) {
			cober$a_indices[[i]][, 2] <- specnumber(
				cober$samples[[i]][, 2:(ncol(cober$samples[[i]]) - 22)]
			)
		}
		# Shannon diversity
		for (i in 1:keyvals$n_levels) {
			cober$a_indices[[i]][, 3] <- diversity(
				cober$samples[[i]][, 2:(ncol(cober$samples[[i]]) - 22)],
				index = "shannon"
			)
		}
		# Simpson diversity
		for (i in 1:keyvals$n_levels) {
			cober$a_indices[[i]][, 4] <- diversity(
				cober$samples[[i]][, 2:(ncol(cober$samples[[i]]) - 22)],
				index = "simpson"
			)
		}
		# inverse Simpson
		for (i in 1:keyvals$n_levels) {
			cober$a_indices[[i]][, 5] <- diversity(
				cober$samples[[i]][, 2:(ncol(cober$samples[[i]]) - 22)],
				index = "invsimpson"
			)
		}
	## carnon
	
	# species richness
	for (i in 1:keyvals$n_levels) {
		carnon$a_indices[[i]][, 2] <- specnumber(
			carnon$samples[[i]][, 2:(ncol(carnon$samples[[i]]) - 22)]
		)
	}
	# Shannon diversity
	for (i in 1:keyvals$n_levels) {
		carnon$a_indices[[i]][, 3] <- diversity(
			carnon$samples[[i]][, 2:(ncol(carnon$samples[[i]]) - 22)],
			index = "shannon"
		)
	}
	# Simpson diversity
	for (i in 1:keyvals$n_levels) {
		carnon$a_indices[[i]][, 4] <- diversity(
			carnon$samples[[i]][, 2:(ncol(carnon$samples[[i]]) - 22)],
			index = "simpson"
		)
	}
	# inverse Simpson
	for (i in 1:keyvals$n_levels) {
		carnon$a_indices[[i]][, 5] <- diversity(
			carnon$samples[[i]][, 2:(ncol(carnon$samples[[i]]) - 22)],
			index = "invsimpson"
		)
	}
	
	## red
	
		# species richness
		for (i in 1:keyvals$n_levels) {
			red$a_indices[[i]][, 2] <- specnumber(
				red$samples[[i]][, 2:(ncol(red$samples[[i]]) - 22)]
			)
		}
		# Shannon diversity
		for (i in 1:keyvals$n_levels) {
			red$a_indices[[i]][, 3] <- diversity(
				red$samples[[i]][, 2:(ncol(red$samples[[i]]) - 22)],
				index = "shannon"
			)
		}
		# Simpson diversity
		for (i in 1:keyvals$n_levels) {
			red$a_indices[[i]][, 4] <- diversity(
				red$samples[[i]][, 2:(ncol(red$samples[[i]]) - 22)],
				index = "simpson"
			)
		}
		# inverse Simpson
		for (i in 1:keyvals$n_levels) {
			red$a_indices[[i]][, 5] <- diversity(
				red$samples[[i]][, 2:(ncol(red$samples[[i]]) - 22)],
				index = "invsimpson"
			)
		}
}