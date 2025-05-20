# Dependencies ----

	{
		# install.packages("vegan")
		library(vegan)
	
		# install.packages("ggplot2")
		library(ggplot2)
	}

# Data Import ----

	{
		# import sample metadata
		bact_points <- list()
		bact_points$all <- read.csv(
			file = "raw_data/sample_metadata_blanks_removed.csv"
		)
		# import ASV data across all taxonomic levels
#		asv_table <- read.csv(
#			"raw_data/asv_table.csv"
#		)
		# import ASV data filtered by each taxonomic level, and bundles them into a list
		bact_data <- list()
		bact_data$ASVs <- list(
			level_2 = read.csv("raw_data/PCR/level_2.csv"),
			level_3 = read.csv("raw_data/PCR/level_3.csv"),
			level_4 = read.csv("raw_data/PCR/level_4.csv"),
			level_5 = read.csv("raw_data/PCR/level_5.csv"),
			level_6 = read.csv("raw_data/PCR/level_6.csv")
		)
		# rename the index column to match "sample_metadata"
		for (i in 1:length(bact_data$ASVs)) {
			colnames(bact_data$ASVs[[i]])[1] <- "SampleID"
		}
		# check output
		print(bact_data$ASVs[[1]][1:4,1:4])
	}

# Data Formatting ----

	# Initialise and order factor columns
	# {
	# 	sample_metadata$SampleID <- as.factor(
	# 		sample_metadata$SampleID
	# 	)
	# 	sample_metadata$Sample <- as.factor(
	# 		
	# 	)


# Parity Checks ----
	
	# Same number of samples in each taxonomic level?
	{
		# initialise parity state vector
		parity <- logical()
		# pairwise comparisons of number of samples in each level
		for (i in 1:length(bact_data$ASVs)) {
			parity[i] <- all.equal(
				nrow(bact_data$ASVs[["level_2"]]),
				nrow(bact_data$ASVs[[i]])
			)
		}
		# result (should all be TRUE)
		print(parity)
		# clean up temporary objects
		rm(parity)
	}
		
	# Same sample IDs at each taxonomic level?
	{	
		# initialise parity state vector
		parity <- logical()
		# pairwise comparisons of sample IDs in each level
		for (i in 1:length(bact_data$ASVs)) {
			parity[i] <- all.equal(
				bact_data$ASVs[["level_2"]]["SampleID"],
				bact_data$ASVs[[i]]["SampleID"]
			)
		}
		# result (should all be TRUE)
		print(parity)
		# clean up temporary objects
		rm(parity)
	}
		
	# Confirmed: Each taxonomic level contains identical sample numbers and names

# Metadata ----
	
	{
		# number of taxonomic levels examined in this project
		keyvals$n_levels <- length(bact_data$ASVs)
		# total numbers of samples in this project
		keyvals$n_samples_tot <- nrow(bact_points$all)
		# number of samples actually sequenced
		keyvals$n_samples_seq <- nrow(bact_data$ASVs[["level_2"]])
	}

	# number of operational taxonomic units at each taxonomic level
	{
		# matrix to store the number of OTUs at each level
		keyvals$n_OTUs <- matrix(
			0, keyvals$n_levels, 1
		)
		# fill "n_OTUs" with the number of OTUs at each level
		for (i in 1:keyvals$n_levels) {
			# subtract 1 to discount the "SampleID" column
			keyvals$n_OTUs[i] <- ncol(bact_data$ASVs[[i]]) - 1
		}
		# label "n_OTUs" rows and cols
		rownames(keyvals$n_OTUs) <- names(bact_data$ASVs)
		colnames(keyvals$n_OTUs) <- "OTU_Number"
	}

# New Lists ----
	## creates a copy of ASVs where sample IDs are stored in the row names, instead of being mixed in with the data
	## separates data by catchment
	{
	# Copy of "ASVs" with no "SampleID" column
		# copy "ASVs", excluding the "SampleID" column
		bact_data$rowlabelled <- list()
		for (i in 1:keyvals$n_levels) {
			bact_data$rowlabelled[[i]] <- bact_data$ASVs[[i]]
			rownames(bact_data$rowlabelled[[i]]) <- bact_data$rowlabelled[[i]][, "SampleID"]
			bact_data$rowlabelled[[i]] <- bact_data$rowlabelled[[i]][, !colnames(bact_data$rowlabelled[[i]]) %in% "SampleID"]
		}
		# copy taxonomic level names
		names(bact_data$rowlabelled) <- names(bact_data$ASVs)
	}
	
	# Sample-Species Database
#	{
#		# initialise list
#		sample_species <- list()
#		for (i in 1:n_levels) {
#			temp_list <- list()
#			# Temporary list to store data frames
#			for (j in 1:n_samples_seq) {
#				# Get species and sample information
#				species <- colnames(
#					ASVs_unlabelled[[i]][
#						which(ASVs_unlabelled[[i]][j, ] != 0)
#					]
#				)
#				sample <- rep(
#					x = ASVs[["level_2"]]$SampleID[j],
#					times = length(species)
#				)
#				# Store as a data frame in the list
#				temp_list[[j]] <- data.frame(
#					species = species,
#					sample = sample
#				)
#			}
#			# Combine all data frames for taxonomic level i
#			sample_species[[i]] <- do.call(rbind, temp_list)
#		}
#		# clean up temporary objects
#		rm(temp_list, sample, species)
#	}
	
	# Sample factor database
	# {
	# 	# initialise list
	# 	test <- list()
	# 		for (j in 1:ncol(sample_metadata)) {
	# 			test[[i]] <-
	# 	for (i in 1:n_levels) {
	# 		}
	# 	}
	# }

# Calculating Alpha-Diversity ----

	# Create Results Storage List
	{
		# template to store results for a single taxonomic level
		template <- data.frame(
			SampleID = bact_data$ASVs[["level_2"]][,"SampleID"],
			Richness = rep(0, keyvals$n_samples_seq),
			Shannon = rep(0, keyvals$n_samples_seq),
			Simpson = rep(0, keyvals$n_samples_seq),
			Invsimpson = rep(0, keyvals$n_samples_seq)
		)
		# copy template for each taxonomic level and bundle into a list
		bact_data$a_indices <- vector(
			mode = "list",
			length = keyvals$n_levels
		)
		bact_data$a_indices <- lapply(
			bact_data$a_indices, function(x) template
		)
		names(bact_data$a_indices) <- names(bact_data$ASVs)
		# clean up temporary objects
		rm(template)
	}

	# Calculate Per-Sample Alpha-Diversity Indices
	{
		# species richness
		for (i in 1:keyvals$n_levels) {
			bact_data$a_indices[[i]][,"Richness"] <- specnumber(
				bact_data$rowlabelled[[i]]
			)
		}
		# Shannon diversity
		for (i in 1:keyvals$n_levels) {
			bact_data$a_indices[[i]][,"Shannon"] <- diversity(
				bact_data$rowlabelled[[i]],
				index = "shannon"
			)
		}
		# Simpson diversity
		for (i in 1:keyvals$n_levels) {
			bact_data$a_indices[[i]][,"Simpson"] <- diversity(
				bact_data$rowlabelled[[i]],
				index = "simpson"
			)
		}
		# inverse Simpson
		for (i in 1:keyvals$n_levels) {
			bact_data$a_indices[[i]][,"Invsimpson"] <- diversity(
				bact_data$rowlabelled[[i]],
				index = "invsimpson"
			)
		}
		# check results
		print(bact_data$a_indices[[1]][1:4,1:5])
		print(bact_data$a_indices[[2]][1:4,1:5])
	}

# Calculating Beta-Diversity ----

	# Generate Bray-Curtis Distance Matrix (locally-generated)
	{
		# initialise storage list
		distances <- list()
		# calculate distance matrices for each level
		for (i in 1:keyvals$n_levels) {
			distances$all[[i]] <- metaMDS(
				comm = bact_data$rowlabelled[[i]],
				distance = "bray",
				k = 2,
				trymax = 100
			)
		}
		# stress values < 0.2?
		for (i in 1:keyvals$n_levels) {
			print(bact_data$distances[[i]]$stress)
		}
	}
	
# Prepare For Plotting ----
	
	# Combine Alpha Metrics with Metadata
	{
		for (i in 1:keyvals$n_levels) {
			merged$alpha.bact_point[[i]] <- merge(
				x = bact_data$a_indices[[i]],
				y = bact_points$all,
				by = "SampleID"
			)
		}
	}
	
	# Combine Bray-Curtis Distance of Each Sample with Metadata
	{
		# Extract and label "points" coordinates
		{
			# transfer coordinates of each sample to storage list
			for (i in 1:keyvals$n_levels) {
				merged$distances.IDs[[i]] <- as.data.frame(
					distances$all[[i]][["points"]]
				)
			}
			# add sample name data
			for (i in 1:keyvals$n_levels) {
				merged$distances.IDs[[i]] <- cbind(
					bact_data$ASVs[["level_2"]][, "SampleID"],
					merged$distances.IDs[[i]]
				)
				# rename label columns to match "sample_metadata"
				colnames(merged$distances.IDs[[i]])[1] <- "SampleID"
			}
			
		}
		
		# Merge sample coordinates and metadata
		{
			for (i in 1:keyvals$n_levels) {
				# merge at each taxonomic level
				merged$distances.metadata[[i]] <- merge(
					x = merged$distances.IDs[[i]],
					y = bact_points$all,
					by = "SampleID"
				)
			}
			# name taxonomic levels
			names(merged$distances.metadata) <- names(bact_data$ASVs)
		}
		# clean up temporary objects
#		rm(merged$distances.IDs)
	}

#	# open external plotting window
#	{
#		  windows(width = 10, height = 8)
#		if (.Platform$OS.type == "windows") {
#		} else if (.Platform$OS.type == "unix") {
#		  quartz(width = 10, height = 8)  # For macOS
#		} else {
#		  x11(width = 10, height = 8)  # For Linux/Unix
#		}
#	}

# Create a Proper Distance Matrix

	{
		dist_proper <- ASVs
		for (i in 1:n_levels) {
			rownames(dist_proper[[i]]) <- dist_proper[[i]]$SampleID
			dist_proper[[i]] <- dist_proper[[i]][, -1]
			dist_proper[[i]] <- as.matrix(dist_proper[[i]])
		}
	}


# Errors / Problems / Questions ----
			
	# fairly high stress
	# "braycurtis.tsv" distance matrix not working with metaMDS? (see section immediately below)
	# need "metadata" referenced in ("analysis_R_script") to run "analysis_R_script"?
	# If nMDS draws points by species, how can sites be differentiated?
		# potentially multiple sites per species - one site cannot own one species
	# How to reliably match point to metadata? Use "points" element?
				
# Next Steps ----
		
	# Colour bar chart by species, separate by river
	# Filter nMDS by catchment
	
# Notepad ----

	# box plots of a-div metrics
	# grouped by factor level
		# distance to mouth
	# qualitativiely state the level of pollution of each river

	# good graphs
		# box whisker by month

# Next Week ----
	
	# Introduction
	# Methods
	# How to integrate EA environmental data with sample env. data?
		# Sample env data, EA metal data
	# Perform statistical tests
		# gradient in EA env. data? (interesting because it may serve as an indicator of pollution, AND may evoke bacterial change)
		# gradients in sample env. data?
	# Main graphs:
		# 
# relative binning of sites into "high, medium, low"
# abundance
	# 
# alpha-diversity
	# box plot colour coded alpha-diversity y-axis
	# work out whether I need to account for seasonality?
# beta-diversity
	# statistical PERMANOVA / ADONIS
	# NDMS plot (separated by catcvhemnt and colour coded by pollution level)
# pH, distance to mouth, metals
# flagship species - which are more abundant at high pollution, are they tolerant to multiple stresses? How to evaluate the strength of metal correlation if a particualr speices is multi-stress tolerant, and are there unexpected speices that may have adapted in situ
# 