# This pipeline:

	# - Imports and pairs bacterial metadata and concentration data
	# - Puts the data to wide format
	# - Filters aggressively to maximise number of complete rows (removing metals from consideration)
	
	# - Runs PCA on the metal data

	# - Preps for PERMANOVA
	# - If there are multiple PCs per Catchment/Site/Month combo, collapse and aggregate to one per season/month
	# - Perform linear PERMANOVA

{	# Import Metal Sampling Point Metadata ----
	
	# - imports metadata from all metal sampling points in the area
	# - filters out any that are not paired to a bacteria point
	
	# create list	
		metal_points <- list()
	# read in raw data as separate lists per river
		for (i in 1:length(keyvals$catchments)) {
			metal_points[[i]] <- read.csv(
				paste0(
					"raw_data/sampling_points/sp_", keyvals$catchments[i], ".csv"
				),
				header = TRUE
			)
		}
	# combine the data into a single list
		metal_points <- do.call(
			rbind, metal_points
		)
	# make sure "Catchment" is a factor
		metal_points$Catchment <- as.factor(
			metal_points$Catchment
		)
	# keep only metal points with a bacterial point twin
		metal_points <- metal_points[
			!is.na(metal_points$eDNA_point)
			,
		]
}

{	# Import and Format Metal Sampling Point Concentration Data ----

	# - imports metal concentration data from all metal sampling points in the area
	# - filters out data from non-paired metal points based on the section above
	# - filters out non-metal data
	# - transposes data to wide format (each sampling point-month combo gets its own row)
	
	# create master list
		metal_data <- list()
	# import raw data
		for (i in keyvals$years[1]:keyvals$years[2]) {
			metal_data[[i - (keyvals$years[1] - 1)]] <- read.csv(
				paste0(
					"raw_data/EA_water_data/DC-", i, ".csv"
				)
			)
		}
	# rbind each list item together
		metal_data <- do.call(
			rbind, metal_data
		)
	# remove results from un-paired points
		metal_data <- subset(
			x = metal_data,
			subset = sample.samplingPoint.notation %in% metal_points$notation
		)
	# tidy up columns
		# create separate year, month, time columns
			metal_data[c('year', 'month', 'time')] <- str_split_fixed(
				metal_data$sample.sampleDateTime, '-', 3
			)
		# create "year_month" column
			metal_data$point_year_month <- paste(
				metal_data$sample.samplingPoint.notation, metal_data$year, metal_data$month,
				sep = "_"
			)
		# create "measure_unit" column
			metal_data$measure_unit <- paste(
				metal_data$determinand.definition, metal_data$determinand.unit.label,
				sep="_"
			)
		# tidy up "measurement units" contents
			metal_data$measure_unit <- metal_data$measure_unit %>%
				gsub("\\s+", "", .) %>%
				gsub("-", "_", .) %>%
				gsub("/", ".", .)
	# filter out non-metal results
		# remove unnecessary results
			metal_data <- subset(
				x = metal_data,
				subset = determinand.definition %in% keyvals$parameters
			)
		# remove all but essential columns
			metal_data <- metal_data %>%
				select(
					result,
					point_year_month,
					measure_unit
				)
	# export
		write_csv(
			metal_data,
			"outputs/tables/metal_data_long.csv"
		)
	# transpose to wide format
		metal_data <- reshape(
			data = metal_data,
			idvar = "point_year_month",
			timevar = "measure_unit",
			direction = "wide"
		)
	# tidy up wide column names
		# separate point_year_month
			metal_data <- cbind(
				as.data.frame(
					str_split_fixed(
						string = metal_data$point_year_month,
						pattern = "_",
						n = 3
					)
				),
				metal_data
			)
		# name new columns
			colnames(metal_data)[1] <- "notation"
			colnames(metal_data)[2] <- "year"
			colnames(metal_data)[3] <- "month"
		# remove now-redundant point_year_month column
			metal_data <- metal_data %>%
				select(
					-point_year_month
				)
		# Clean up metal column names by overwriting with names from keyvals
			for (i in keyvals$metals) {
				# Find column index where the name appears (partial match)
					match_index <- grep(
						i,
						colnames(metal_data),
						ignore.case = TRUE
					)
				# If a match is found, rename the column
					if (length(match_index) == 1) {
						colnames(metal_data)[match_index] <- i
					}
			}
	# add sampling point metadata from metal_points
		metal_wide <- merge(
			x = metal_points %>%
				select(
					Catchment,
					notation,
					eDNA_point
				),
			y = metal_data,
			by = "notation"
		)
	# convert "eDNA_point" to "Site" to match bact_data
		metal_wide <- metal_wide %>% rename(Site = eDNA_point)
	# remove intermediate dataframes
		rm(
			metal_data,
			metal_points
		)
	# add "season" factor column
		metal_wide <- metal_wide %>% 
			mutate(
				Season = case_when(
					as.numeric(month) %in% c(3, 4, 5)   ~ "Spring",
					as.numeric(month) %in% c(6, 7, 8)   ~ "Summer",
					as.numeric(month) %in% c(9, 10, 11) ~ "Autumn",
					as.numeric(month) %in% c(12, 1, 2)  ~ "Winter"
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
				.after = 4
			)
	# Export wide-format metal data as csv file
		write_csv(
			x = metal_wide,
			file = "outputs/tables/metal_data_wide.csv"
		)
}

{ # Shapiro-Wilk test on all raw metal values ----
	# Justifies the use of medians
		norm_tests <- list()
		norm_tests$rawmetal <- lapply(keyvals$metals, function(m) {
			x <- na.omit(metal_wide[[m]])
			n <- length(x)
			
			if (n < 3 || var(x) == 0) {
				data.frame(Metal = m, W = NA, p.value = NA, n = n,
							  Note = if (n < 3) "<3 obs" else "Constant values")
			} else {
				if (n > 5000) x <- sample(x, 5000)        # optional
				s <- shapiro.test(x)
				data.frame(Metal = m, W = unname(s$statistic),
							  p.value = s$p.value, n = n, Note = "")
			}
		}) |>
			dplyr::bind_rows()
	# export table
		write.csv(
			x = norm_tests$rawmetal,
			file = "outputs/tables/distest_rawmetal.xlsx"
		)
}

{ # Summarise metal data per site/season (median concs, sample counts), and overall (sample count vs NA count) ----
	metal_medians <- metal_wide %>% 
		group_by(Catchment, notation, Site, Season) %>% 
	 	summarise(
	 		across(
	 			all_of(keyvals$metals),
	 			~median(.x, na.rm = TRUE),
	 			.names = "{.col}"
	 		),
	 		.groups = "drop"
	 	)
	 
	 metal_iqr <- metal_wide %>%
	 	group_by(Catchment, notation, Site, Season) %>%
	 	summarise(
	 		across(
	 			all_of(keyvals$metals),
	 			~IQR(.x, na.rm = TRUE),
	 			.names = "{.col}"
	 		),
	 		.groups = "drop"
	 	)
	 
	 metal_counts <- metal_wide %>% 
	 	group_by(Catchment, notation, Site, Season) %>% 
	 	summarise(
	 		across(
	 			all_of(keyvals$metals),
	 			~sum(!is.na(.x)),
	 			.names = "{.col}"
	 		),
	 		.groups = "drop"
	 	)
	 	
metal_summary <- metal_wide %>%
  pivot_longer(
    cols = all_of(keyvals$metals),
    names_to = "Metal",
    values_to = "Value"
  ) %>%
  group_by(Metal) %>%
  summarise(
    Actual = sum(!is.na(Value)),
    NA_count = sum(is.na(Value)),
    Percentage = mean(!is.na(Value)) * 100
  )

complete_percentage <- metal_wide %>%
  mutate(is_complete = complete.cases(select(., all_of(keyvals$metals)))) %>%
  summarise(percent_complete = mean(is_complete) * 100) %>%
  pull()

}

{ # filter poorly-sampled metals to maximise the number of usable complete records ----
	# remove metals with far more NA values than actual values across all metal samples
	# increases number of complete measurements from 177 to 235
		metal_wide <- metal_wide %>%
			select(
				-Aluminium,
				-Barium,
				-Lithium
			)
		metal_medians <- metal_medians %>%
			select(
				-Aluminium,
				-Barium,
				-Lithium
			)
		metal_iqr <- metal_iqr %>%
			select(
				-Aluminium,
				-Barium,
				-Lithium
			)
		metal_counts <- metal_counts %>%
			select(
				-Aluminium,
				-Barium,
				-Lithium
			)
	# remove metals that don't vary much
	# number of complete measurements stays the same
		metal_wide <- metal_wide %>%
			select(
				-Chromium,
				-Boron,
				-Lead
			)
		metal_medians <- metal_medians %>%
			select(
				-Chromium,
				-Boron,
				-Lead
			)
		metal_iqr <- metal_iqr %>%
			select(
				-Chromium,
				-Boron,
				-Lead
			)
		metal_counts <- metal_counts %>%
			select(
				-Chromium,
				-Boron,
				-Lead
			)
	# Remove Arsenic to allow more Hayle points to be included
	# increases complete measurements from 235 to 341
		metal_wide <- metal_wide %>%
			select(
				-Arsenic
			)
		metal_medians <- metal_medians %>%
			select(
				-Arsenic
			)
		metal_iqr <- metal_iqr %>%
			select(
				-Arsenic
			)
		metal_counts <- metal_counts %>%
			select(
				-Arsenic
			)	
}

{ # Remove rows with NA values (and zero values in metal_counts) ----
		metal_wide <- metal_wide %>% drop_na
		metal_medians <- metal_medians %>% drop_na()
		metal_iqr <- metal_iqr %>% drop_na()
		metal_counts <- metal_counts %>% 
			filter(if_all(everything(), ~ .x != 0))
}

{ # cor test the retained metals ----
metal_wide %>%
  select("Cadmium", "Iron", "Nickel", "Copper", "Manganese", "Zinc", "Magnesium") %>%
  cor(method = "spearman", use = "pairwise.complete.obs") %>%
		as.data.frame() %>%
		write_csv(., "outputs/tables/corr_table_final.csv")
}

{ # Export the dataframes ----
		write_csv(
			x = metal_medians,
			file = "outputs/tables/metal_medians.csv"
		)
		write_csv(
			x = metal_iqr,
			file = "outputs/tables/metal_iqr.csv"
		)
		write_csv(
			x = metal_counts,
			file = "outputs/tables/metal_counts.csv"
		)
		write_csv(
			x = metal_wide,
			file = "outputs/tables/metal_wide_trimmed.csv"
		)
}

{ # Run naive PCA and PERMANOVA on metal data (significant over-dispersion) ----
	# factorise the data
		meta <- metal_wide %>% 
			mutate(
				Catchment = factor(Catchment),
				Site      = factor(Site),
				Season    = factor(Season,
										 levels = c("Spring","Summer","Autumn","Winter"))
			)
	# Isolate the metal columns
		metal_cols <- c("Cadmium","Iron","Nickel",
							 "Copper","Manganese","Zinc","Magnesium")
		metals <- meta %>% select(all_of(metal_cols))
	# Standardise metals to remove unit size effects
		metals_z <- scale(metals)               # subtract mean, divide by sd
	# Run Naive PCA on metals
		pca_res <- prcomp(metals_z, center = FALSE, scale = FALSE)
	# Cumulative variance explained by the PCs
		summary(pca_res)$importance[2:3,1:5]
	# Choose how many PCs to keep
		keep <- 3
	# Pull PC scores for those axes
		scores <- as.data.frame(pca_res$x[, 1:keep]) %>% 
			setNames(paste0("PC", 1:keep))
	# Add scores to the metadata table
		meta_pc <- bind_cols(meta, scores)
	# Aggregate > 1 measurement per Catchment/Site/Year/Season
		meta_pc <- meta_pc %>% 
			group_by(Catchment, Site, year, Season) %>%
			summarise(across(starts_with("PC"), median),
						 .groups = "drop")
	# Check loading values of each PC
		naive_loadings <- pca_res$rotation[, 1:keep]  # weights of each metal
		print(round(naive_loadings, 2))               # keep two decimals
	# Naive PERMANOVA - whether the metal mixes cluster between factor levels
		perm <- adonis2(
			meta_pc[, paste0("PC", 1:keep)] ~
				Catchment + Site %in% Catchment + year %in% Catchment + Season,
			data = meta_pc,
			by = "terms",
			method = "euclidean",
			strata = meta_pc$Site
		)
		print(perm)
	# export
		write_csv(
			meta_pc,
			"outputs/tables/meta_pc.csv"
		)
		write_csv(
			meta,
			"outputs/tables/meta.csv"
		)
	# Convert matrix to dataframe with row names as a column
		naive_loadings_df <- as.data.frame(naive_loadings) %>% 
			rownames_to_column(var = "Metal")  # Convert row names to a column
	# Write to CSV
		write_csv(
			naive_loadings_df,
			"outputs/tables/naive_loadings_df.csv"
		)
	# convert to long format for comparison against iterative balanced PCA in the next step
		naive_long <- naive_loadings_df %>%
			pivot_longer(
				cols = -Metal,
				names_to = "PC",
				values_to = "orig_loading"
			)
}

{ # test naive PC values for over-dispersion across each factor ----
	metal_dist <- dist(meta_pc %>% select(PC1:paste0("PC", keep)))
	# over-dispersion by catchment
		permutest(
			betadisper(
				metal_dist,
				group = meta_pc$Catchment
			),
			permutations = 999
		)
	# over-dispersion by site in catchment
		permutest(
			betadisper(
				metal_dist,
				group = meta_pc$Site
			),
			permutations = how(
				nperm = 999,
				within = Within(type = "free"),
				plots = Plots(strata = meta_pc$Catchment)
			)
		)
	# over-dispersion by year in catchment
		permutest(
			betadisper(
				metal_dist,
				group = meta_pc$year
			),
			permutations = how(
				nperm = 999,
				within = Within(type = "free"),
				plots = Plots(strata = meta_pc$Catchment)
			)
		)
	# over-dispersion by season
		permutest(
			betadisper(
				metal_dist,
				group = meta_pc$Season
			),
			permutations = 999
		)
}

{ # post-hoc to detect exactly HOW metal scores vary with each factor ----
		pcs <- meta_pc %>%
			select(PC1:paste0("PC", keep)) %>%
			names()
	
	# ── fixed helper
	kw_boxplot <- function(data, pc, factor){
		
		# 1) run Kruskal + Dunn grouping
		kw_out <- kruskal(
			y   = data[[pc]],
			trt = as.factor(data[[factor]]),
			group     = TRUE,
			p.adj     = "bonferroni"
		)
		
		revo <- kruskal(
			y   = data[[pc]],
			trt = as.factor(data[[factor]]),
			group     = FALSE,
			p.adj     = "bonferroni"
		)
		
		# ---- A.  append key stats to kw_log
		cat(pc, " ", factor, "\n")
		print(revo$statistics)
		print(revo$comparison)
		
		grp_df <- kw_out$groups %>%
			rownames_to_column(var = factor) %>%
			rename(letter = groups)
		
		plot_df <- data %>%
			mutate(across(all_of(factor), as.factor)) %>%
			left_join(grp_df, by = factor)
		
		pal <- scales::hue_pal()(n_distinct(plot_df$letter))
		names(pal) <- sort(unique(plot_df$letter))
		
		pal_fill   <- scales::alpha(pal, 0.10)   # 40 % opaque
		names(pal_fill) <- names(pal)            # keep the same letter names
		
		p <- ggplot(plot_df, aes_string(x = factor, y = pc, fill = "letter", colour = "letter")) +
			geom_boxplot(alpha = 0.1, show.legend = FALSE, outlier.shape = NA, size = 0.6) +   # outlines coloured, fill semi-transparent
			geom_jitter(width = 0.15, alpha = 0.7, size = 1, shape = 1) +
#			scale_x_discrete(limits = levels(meta_pc[[factor]])) +
			scale_fill_manual(values = pal_fill) +
			scale_colour_manual(values = pal) +              # fully opaque for borders & text
			geom_text(
				data = grp_df,
				aes_string(x = factor,
							  y = max(plot_df[[pc]], na.rm = TRUE) * 1.3,
							  label = "letter",
							  colour = "letter"),
				show.legend = FALSE
			) +
			theme_classic() +
			theme(
				plot.subtitle = element_text(size = 9),
				legend.position = "none"
			)
		
		print(p)
		
		ggsave(
			filename = paste0("outputs/plots/final/", pc, "_by_", factor, ".png"),  # Dynamic filename
			plot = p,
			width = 3,
			height = 5,
			units = "in",
			dpi = 1000
		)
	}
}

{ # post-hoc graphs for all PCs and factors ----
	
	# Each standalone factor across all PC axes
	walk(pcs, ~ print(kw_boxplot(meta_pc, pc = .x, factor = "Catchment")))
	walk(pcs, ~ print(kw_boxplot(meta_pc, pc = .x, factor = "Season")))
	
	# sites within a specific catchment
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Carnon"),
			pc     = pc,
			factor = "Site")
		print(p)
	})
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Cober"),
			pc     = pc,
			factor = "Site")
		print(p)
	})
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Hayle"),
			pc     = pc,
			factor = "Site")
		print(p)
	})
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Red"),
			pc     = pc,
			factor = "Site")
		print(p)
	})
	
	# years within a specific catchment
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Carnon"),
			pc     = pc,
			factor = "year")
		print(p)
	})
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Cober"),
			pc     = pc,
			factor = "year")
		print(p)
	})
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Hayle"),
			pc     = pc,
			factor = "year")
		print(p)
	})
	walk(pcs, function(pc){
		p <- kw_boxplot(
			data   = meta_pc %>% filter(Catchment == "Red"),
			pc     = pc,
			factor = "year")
		print(p)
	})
}

{ # Catchment-level nMDS plot ----
	## pick easily recognisable symbols
	shape_vals <- c(
		"Carnon" = 16,
		"Red"    = 17,
		"Cober"  = 15,
		"Hayle"  = 18
	)
	
	# ----- 2. build the distance matrix on the four PCs
	# PCs are already scaled, so Euclidean distance is fine
	dist_mat <- dist(
		meta_pc %>% select(PC1:paste0("PC", keep)),
		method = "euclidean")
	
	# ----- 3. run non-metric MDS (two dimensions)
	set.seed(123)                               # gives reproducible solution
	fit <- metaMDS(dist_mat, k = 2, trymax = 100)
	
	# ----- 4. combine ordination scores with the metadata
	plot_df <- meta_pc %>%                     # 341 rows
		mutate(MDS1 = fit$points[, 1],
				 MDS2 = fit$points[, 2])
	
	plot_df$Catchment <- factor(plot_df$Catchment, levels = names(shape_vals))
	
	# ----- 5. draw the map
	p <- ggplot(plot_df, aes(MDS1, MDS2,
									 colour = Catchment,
									 shape = Catchment)) +
		geom_point(size = 3, alpha = 0.8) +
		scale_shape_manual(values = shape_vals) +
		theme_classic() +
		labs(title = NULL,
			  x = "nMDS Axis 1", y = "nMDS Axis 2")
	p <- p + ggforce::geom_mark_hull(aes(fill = Catchment),
										 concavity = 5, linetype = 0,
										 alpha = 0.08)
	ggsave(
		filename = paste0("outputs/plots/final/Catchments.png"),  # Dynamic filename
		plot = p,
		width = 9,
		height = 5,
		units = "in",
		dpi = 1000
	)
	p
}

{ # Aggregate PC values by season, collapse different years together ----
	meta_pc_agg <- meta_pc %>%
		group_by(
			Site,
			Season
		) %>%
		summarise(
			across(
				starts_with("PC"),
				\(x) median(x, na.rm = TRUE)
			),
			.groups = "drop"
		)
}
