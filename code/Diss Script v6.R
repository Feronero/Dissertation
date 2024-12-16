# Dependencies ----

	{
		# set working directory
		setwd(
			"~/Documents/University/Year 3/Term 1/Dissertation/R"
		)
		# install.packages("vegan")
		library(vegan)
	
		# install.packages("ggplot2")
		library(ggplot2)
	}

# Data Import ----

	{
		# import sample metadata
		sample_metadata <- read.csv(
			file = "Data/sample_metadata_ordered.csv"
		)
		# import ASV data across all taxonomic levels
		asv_table <- read.csv(
			"Data/asv_table.csv"
		)
		# import ASV data filtered by each taxonomic level, and bundles them into a list
		ASVs <- list(
			level_2 = read.csv("Data/level_2.csv"),
			level_3 = read.csv("Data/level_3.csv"),
			level_4 = read.csv("Data/level_4.csv"),
			level_5 = read.csv("Data/level_5.csv"),
			level_6 = read.csv("Data/level_6.csv")
		)
		# rename the index column to match "sample_metadata"
		for (i in 1:length(ASVs)) {
			colnames(ASVs[[i]])[1] <- "SampleID"
		}
		# check output
		print(ASVs[[1]][1:4,1:4])
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
		for (i in 1:length(ASVs)) {
			parity[i] <- all.equal(
				nrow(ASVs[["level_2"]]),
				nrow(ASVs[[i]])
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
		for (i in 1:length(ASVs)) {
			parity[i] <- all.equal(
				ASVs[["level_2"]]["SampleID"],
				ASVs[[i]]["SampleID"]
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
		n_levels <- length(ASVs)
		# total numbers of samples in this project
		n_samples_tot <- nrow(sample_metadata)
		# number of samples actually sequenced
		n_samples_seq <- nrow(ASVs[["level_2"]])
	}

	# number of operational taxonomic units at each taxonomic level
	{
		# matrix to store the number of OTUs at each level
		n_OTUs <- matrix(
			0, n_levels, 1
		)
		# fill "n_OTUs" with the number of OTUs at each level
		for (i in 1:n_levels) {
			# subtract 1 to discount the "SampleID" column
			n_OTUs[i] <- ncol(ASVs[[i]]) - 1
		}
		# label "n_OTUs" rows and cols
		rownames(n_OTUs) <- names(ASVs)
		colnames(n_OTUs) <- "OTU_Number"
	}

# New Lists ----

	# Copy of "ASVs" with no "SampleID" column
	{
		# initialise list
		ASVs_unlabelled <- vector(
			mode = "list",
			length = n_levels
		)
		# copy "ASVs", excluding the "SampleID" column
		for (i in 1:n_levels) {
			ASVs_unlabelled[[i]] <- ASVs[[i]][, names(ASVs[[i]]) != "SampleID"]
		}
		# copy taxonomic level names
		names(ASVs_unlabelled) <- names(ASVs)
	}
	
	# Sample-Species Database
	{
		# initialise list
		sample_species <- list()
		for (i in 1:n_levels) {
			# Temporary list to store data frames
			temp_list <- list()
			for (j in 1:n_samples_seq) {
				# Get species and sample information
				species <- colnames(
					ASVs_unlabelled[[i]][
						which(ASVs_unlabelled[[i]][j, ] != 0)
					]
				)
				sample <- rep(
					x = ASVs[["level_2"]]$SampleID[j],
					times = length(species)
				)
				# Store as a data frame in the list
				temp_list[[j]] <- data.frame(
					species = species,
					sample = sample
				)
			}
			# Combine all data frames for taxonomic level i
			sample_species[[i]] <- do.call(rbind, temp_list)
		}
		# clean up temporary objects
		rm(temp_list, sample, species)
	}
	
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
			SampleID = ASVs[["level_2"]][,"SampleID"],
			Richness = rep(0, n_samples_seq),
			Shannon = rep(0, n_samples_seq),
			Simpson = rep(0, n_samples_seq),
			Invsimpson = rep(0, n_samples_seq)
		)
		# copy template for each taxonomic level and bundle into a list
		a_indices <- vector(
			mode = "list",
			length = n_levels
		)
		a_indices <- lapply(
			a_indices, function(x) template
		)
		names(a_indices) <- names(ASVs)
		# clean up temporary objects
		rm(template)
	}

	# Calculate Per-Sample Alpha-Diversity Indices
	{
		# species richness
		for (i in 1:n_levels) {
			a_indices[[i]][,"Richness"] <- specnumber(
				ASVs_unlabelled[[i]]
			)
		}
		# Shannon diversity
		for (i in 1:n_levels) {
			a_indices[[i]][,"Shannon"] <- diversity(
				ASVs_unlabelled[[i]],
				index = "shannon"
			)
		}
		# Simpson diversity
		for (i in 1:n_levels) {
			a_indices[[i]][,"Simpson"] <- diversity(
				ASVs_unlabelled[[i]],
				index = "simpson"
			)
		}
		# inverse Simpson
		for (i in 1:n_levels) {
			a_indices[[i]][,"Invsimpson"] <- diversity(
				ASVs_unlabelled[[i]],
				index = "invsimpson"
			)
		}
		# check results
		print(a_indices[[1]][1:4,1:5])
		print(a_indices[[2]][1:4,1:5])
	}

# Calculating Beta-Diversity ----

	# Generate Bray-Curtis Distance Matrix (locally-generated)
	{
		# initialise storage list
		distances <- vector(
			mode = "list",
			length = n_levels
		)
		# name levels
		names(distances) <- names(ASVs)
		# calculate distance matrices for each level
		for (i in 1:n_levels) {
			distances[[i]] <- metaMDS(
				comm = ASVs_unlabelled[[i]],
				distance = "bray",
				k = 2,
				trymax = 100
			)
		}
		# stress values < 0.2?
		for (i in 1:n_levels) {
			print(distances[[i]]$stress)
		}
	}
	

# Prepare For Plotting ----
	
	# Combine Alpha Metrics with Metadata
	{
		coda_alpha <- list()
		for (i in 1:n_levels) {
			coda_alpha[[i]] <- merge(
				x = a_indices[[i]],
				y = sample_metadata,
				by = "SampleID"
			)
		}
	}
	
	# Combine Bray-Curtis Distance of Each Sample with Metadata
	{
		# Extract and label "points" coordinates
		{
			# initialise storage list
			coda <- vector(
				mode = "list",
				length = n_levels
			)
			# transfer coordinates of each sample to storage list
			for (i in 1:n_levels) {
				coda[[i]] <- as.data.frame(
					distances[[i]][["points"]]
				)
			}
			# add sample name data
			for (i in 1:n_levels) {
				coda[[i]] <- cbind(
					ASVs[["level_2"]][, "SampleID"],
					coda[[i]]
				)
				# rename label columns to match "sample_metadata"
				colnames(coda[[i]])[1] <- "SampleID"
			}
		}
		
		# Merge sample coordinates and metadata
		{
			# initialise final storage list
			coda_points <- vector(
				mode = "list",
				length = n_levels
			)
			for (i in 1:n_levels) {
				# merge at each taxonomic level
				coda_points[[i]] <- merge(
					x = coda[[i]],
					y = sample_metadata,
					by = "SampleID"
				)
			}
			# name taxonomic levels
			names(coda_points) <- names(ASVs)
		}
		# clean up temporary objects
		rm(coda)
	}

	# open external plotting window
	{
		if (.Platform$OS.type == "windows") {
		  windows(width = 10, height = 8)
		} else if (.Platform$OS.type == "unix") {
		  quartz(width = 10, height = 8)  # For macOS
		} else {
		  x11(width = 10, height = 8)  # For Linux/Unix
		}
	}
	
# Alpha-Diversity Plots ----

	# box and whisker plots
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

	# box and whisker plots TEMPORARY FIX
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

# My Stacked Bar Chart ----
	
	# input factor
	factor <- "Month"
	# initialise list
	test <- list()
	# encode grouping vector based on input
	for (i in 1:n_levels) {
		# transforms result into factor levels before saving
		grouping_vector <- as.factor(
			
	
	# initialise merge list
	coda_stacked <- vector(
		mode = "list",
		length = n_levels
	)
	for (i in 1:n_levels) {
		# merge at each taxonomic level
		coda_stacked[[i]] <- merge(
			x = ASVs[[i]],
			y = sample_metadata,
			by = "SampleID"
		)
			# add richness column to each merged data frame
			coda_stacked[[i]]$richness <- rep(
				x = 0,
				times = n_samples_seq
			)
	}
	# fill richness column
	for (i in 1:n_levels) {
		for (j in 1:n_samples_seq) (
			coda_stacked[[i]]$richness[j] <- sum(
				# counts non-zero values from "ASV_unlabelled"
				ASVs_unlabelled[[i]][j, ] != 0
			)
		)
	}
	
	for (i in 1:n_levels) {
		ggplot(
			aes(
				x = 1,
				y = 1,
				fill = 1,
			),
			data = 
		)
	}

# Tamsyn Stacked Bar Chart ----

	# selecting the most abundant 20 ASVs
	Bar_OTUs <- read.csv(
		file = "asv_table.csv",
		row.names = 1,
		check.names = F
	)
	Bar_OTUs <- cbind(
		Bar_OTUs, rowSums(Bar_OTUs)
	)
	Bar_OTUs <- Bar_OTUs[order(rowSums(-Bar_OTUs)),]
	Bar_OTUs[,c("rowSums(Bar_OTUs)")] <- list(NULL)
	Bar20 <- head(Bar_OTUs, 20)
	
	tBar20 <- t(Bar20) #transposing data
	mtop20 = melt(tBar20, id.vars = c("ASV"))
	
	mtop20$Population <- metadata$Population
	mtop20$Treatment <- metadata$treatment
	mtop20$Group <- metadata$group
	mtop20$Rep <- metadata$rep
	
	colours = c("#8B0000",  "#008000", "#FF0000","#90EE90",  "#008080" , "#F6AE2D", "#DC134C", "#7B68EE", "#DDA0DD", "#F9ECCC", "#000000",  "#A54657",  "#582630", "#F7EE7F", "#191970", "#679289",  "#33658A","#F1A66A", "#86BBD8", "#4DAA57")
	
	ggplot(mtop20, aes(x = Rep, fill = Var2, y = value)) + 
	  geom_bar(stat = "identity", colour = "black")+ facet_grid(rows= vars(Population), cols=vars(Treatment), scales = "free")+
	  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5), 
	        axis.title.y = element_text(size = 10), legend.title = element_text(size = 12), 
	        legend.text = element_text(size = 10), 
	        axis.text.y = element_text(size = 10)) + 
	  labs(x = "", y = "Abundance (counts)", fill = "OTU")+
	  scale_fill_manual(values = colours)
	
	ggplot(mtop20, aes(x = Rep, y = Var2)) + 
	  geom_point(aes(size= value, fill=Var2), alpha=0.75, shape=21,)+
	  facet_grid(cols=vars(Treatment), rows=vars(Population), scales = "free", space= 'free')+
	  scale_size_continuous(limits = c(1, 9000), range = c(1,10), breaks = c(20,100,500,2500)) +
	  theme(legend.position = "none")+
	  theme(axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5, ), 
	        axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), 
	        legend.text = element_text(size = 8), 
	        axis.text.y = element_text(size = 8)) + 
	  labs(x = "", y = "", fill = "", size= "Counts")
	
	  scale_fill_manual(values = palette30)

# Default nMDS and Stress Plot ----
		
	{
		# stress plots by level
		for (i in 1:n_levels) {
			stressplot(
				object = distances[[i]],
				xlim = range(0:1),
				ylim = range(0:4)
			)
		}
		# nMDS plots by level
		for (i in 1:n_levels) {
			plot(
				x = distances[[i]],
				type = "p",
				main = "NMDS Plot",
				xlim = c(-2.5,2.5),
				ylim = c(-2.5,2.5)
			)
		}
	}
	
# Custom nMDS Plots ----
		
	# plain plot
	for (i in 1:n_levels) {
		 plot(
		 	formula = MDS2 ~ MDS1,
		 	xlim = range(coda_points[[n_levels]][, "MDS1"]),
			ylim = range(coda_points[[n_levels]][, "MDS2"]),
			type = "p",
			data = coda_points[[i]]
		)
	}
	  
	# Colour By Factor
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

# Stacked Bar Chart Trial (with internal data) ----
		
	# Assumptions
		
		# ASV data is interpreted as presence / absence
		# Equal sampling effort between rivers (so bar height can be meaningfully compared)

	# What do I want to show?
		
		# separate bar per river
		# bar height = abundance in arbitrary units
		# bar colour = species
		
		# plot for entire river
		# one bar for each individual sample
		# visually slightly cluster replicate 
		
	# my attempt
	colours = c(
		"#8B0000","#008000","#FF0000","#90EE90",
		"#008080","#F6AE2D","#DC134C","#7B68EE",
		"#DDA0DD","#F9ECCC","#000000","#A54657",
		"#582630","#F7EE7F","#191970","#679289",
		"#33658A","#F1A66A","#86BBD8","#4DAA57"
	),
	ggplot(
		data = mtop20,
		mapping = aes(
			x = Rep,
			fill = Var2,
			y = value
		)
	) + 
	geom_bar(
		stat = "identity",
		colour = "black"
	) +
	facet_grid(
		rows = vars(Population),
		cols = vars(Treatment),
		scales = "free"
	) +
	theme(
		axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5), 
	     axis.title.y = element_text(size = 10),
		legend.title = element_text(size = 12),
		legend.text = element_text(size = 10), 
		axis.text.y = element_text(size = 10)
	) +
	labs(
		x = "",
		y = "Abundance (counts)",
		fill = "OTU"
	) +
	scale_fill_manual(values = colours)
		
# Simple example stacked bar chart ----
library(ggplot2)

# Create example data
data <- data.frame(
  Category = c("A", "A", "B", "B"),
  Subgroup = c("X", "Y", "X", "Y"),
  Value = c(3, 5, 4, 7)
)

# Create a stacked bar chart
ggplot(data, aes(x = Category, y = Value, fill = Subgroup)) +
  geom_bar(stat = "identity") +
  labs(title = "Stacked Bar Chart", x = "Category", y = "Value") +
  theme_minimal()


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

	# Integrate metal data
	
# Notepad ----

	# box plots of a-div metrics
	# grouped by factor level
		# distance to mouth
	# qualitativiely state the level of pollution of each river

	# good graphs
		# box whisker by month