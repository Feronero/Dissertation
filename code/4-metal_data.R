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
	
#	# Create list of each metal and its measure unit for reference
#		metal_data_ref <- data.frame(
#			metals_and_units = unique(
#				metal_data$measure_unit
#			)
#		)
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
	# remove now-unneeded metal_data
		rm(
			metal_data
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
	# # remove "month" column
	# 	metal_wide <- metal_wide %>%
	# 		select(
	# 			-month
	# 		)
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
		write_xlsx(
			x = norm_tests$rawmetal,
			path = "outputs/tables/distest_rawmetal.xlsx"
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
	 	summarise(
	 		across(
	 			all_of(keyvals$metals),
	 			list(
	 				Actual   = ~sum(!is.na(.x)),
	 				NA_count = ~sum(is.na(.x))
	 			),
	 			.names = "{.col}_{.fn}"
	 		)
	 	) %>% 
	 	pivot_longer(
	 		everything(),
	 		names_to   = c("Metal", ".value"),
	 		names_sep  = "_"
	 	)
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
	# Principal Component Analysis (PCA)
	pca_res <- prcomp(metals_z, center = FALSE, scale. = FALSE)
	# Cumulative variance explained by the PCs
	summary(pca_res)$importance[2:3,1:5]
	# Choose how many PCs to keep
	keep <- 4
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
	# test for overdispersion
	d <- dist(meta_pc %>% select(PC1:PC4))
	bd <- betadisper(d, group = meta_pc$Catchment)
	permutest(bd, permutations = 999)
	# 8. Type I PERMANOVA
	perm <- adonis2(
		meta_pc[, paste0("PC", 1:keep)] ~
			Catchment + Site %in% Catchment + year + Season,
		data = meta_pc,
		permutations = 999,
		by = "terms",
		method = "euclidean",
	)
	print(perm)
	# # Shows samples on the two PC axes, coloured by Catchment
	# fviz_pca_ind(pca_res,
	# 				 geom.ind  = "point",
	# 				 pointshape = 21,
	# 				 pointsize  = 2,
	# 				 col.ind   = meta_pc$Catchment) +
	# 	theme_classic()
	# export data
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
}

{ # balanced PERMANOVA sub-sample loop - removes over-dispersion issues ----
	## the smallest river has min_n rows
	min_n <- meta_pc %>% count(Catchment) %>% pull(n) %>% min()
	
	## function that does one balanced draw + PERMANOVA
	one_run <- function(n_sub = min_n) {
		
		# 1. draw n_sub rows from every Catchment
		samp <- meta_pc %>%
			group_by(Catchment) %>%
			slice_sample(n = n_sub) %>%
			ungroup()
		
		# 2. build Euclidean distance on the four PCs
		dist_mat <- dist(samp %>% select(PC1:PC4))
		
		# 3. PERMANOVA with the full formula
		mod <- adonis2(
			dist_mat ~ Catchment + Site %in% Catchment + year + Season,
			data = samp,
			by = "terms",
			permutations = 999
		)
		
		# 4. keep just the rows you care about
		tibble(term = rownames(mod)[1:4],
				 R2   = mod$R2[1:4],
				 p    = mod$`Pr(>F)`[1:4])
	}
	
	## 999 balanced runs
	set.seed(42)
	out <- map_dfr(1:999, ~one_run())
	
	## summarise
	PERM_loop_summary <- out %>%
		group_by(term) %>%
		summarise(median_R2 = median(R2),
					 lo_R2     = quantile(R2, 0.025),
					 hi_R2     = quantile(R2, 0.975),
					 pct_p_lt_0.05 = mean(p < 0.05) * 100,
					 .groups = "drop")
	
	print(PERM_loop_summary)
}

{ # Balanced PCA sub-sample loop - checks whether naive PCA loadings are correct ----
	# reference rotation – enforce metal_cols order
	ref_rot <- naive_loadings_df %>%
		arrange(match(Metal, metal_cols)) %>%   # <<< new line
		column_to_rownames("Metal") %>%
		as.matrix()
	
	# helper unchanged
	align_signs <- function(new_rot, ref_rot) {
		ref      <- ref_rot[rownames(new_rot), , drop = FALSE]
		cor_sign <- sign(diag(cor(new_rot, ref)))
		sweep(new_rot, 2, cor_sign, `*`)
	}
	
	# one balanced bootstrap draw
	one_draw <- function(seed) {
		set.seed(seed)
		
		samp <- meta |>
			group_by(Catchment) |>
			slice_sample(n = min_n) |>
			ungroup()
		
		pca <- prcomp(samp[, metal_cols], center = TRUE, scale. = TRUE)
		
		rot_ok <- align_signs(pca$rotation[, 1:4], ref_rot)
		
		as_tibble(rot_ok, rownames = "Metal") |>
			pivot_longer(-Metal, names_to = "PC", values_to = "loading") |>
			mutate(iter = seed)
	}
	
	# 999 balanced PCAs
	set.seed(42)
	boot_load <- map_dfr(1:999, one_draw)
	
	balanced_load_summary <- boot_load %>%
		group_by(Metal, PC) %>%
		summarise(median_loading = median(loading),
					 lo_99 = quantile(loading, .005),
					 hi_99 = quantile(loading, .995),
					 .groups = "drop")
	
	# *** no second sign_fix block ***
	compare <- balanced_load_summary %>%
		left_join(naive_long, by = c("Metal","PC")) %>%
		mutate(within_99 = orig_loading >= pmin(lo_99, hi_99) &
				 	orig_loading <= pmax(lo_99, hi_99),
				 abs_diff = abs(orig_loading - median_loading))
	
	agreement <- compare %>%
		group_by(PC) %>%
		summarise(correlation = cor(median_loading, orig_loading),
					 rmse        = sqrt(mean((orig_loading - median_loading)^2)),
					 pct_within  = mean(within_99) * 100,
					 .groups = "drop")
	
	write_csv(compare, "outputs/tables/loading_comparison_detail.csv")
	write_csv(agreement,  "outputs/tables/agreement.csv")
}

{ # post-hoc to detect exactly HOW metal scores vary with each factor ----
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
		
		p <- ggplot(plot_df, aes_string(x = factor, y = pc, fill = "letter")) +
			geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
			geom_jitter(width = 0.15, alpha = 0.3, size = 1) +
			scale_fill_manual(values = pal) +
			geom_text(
				data = grp_df,
				aes_string(x = factor,
							  y = max(plot_df[[pc]], na.rm = TRUE) * 1.3,  # Increased from 1.05
							  label = "letter",
							  colour = "letter"),
				show.legend = FALSE
			) +
			scale_colour_manual(values = pal) +
#			labs(title = paste(pc, "by", factor),
#				  subtitle = "Levels sharing a colour/letter are not significantly different\n(Kruskal-Wallis + Bonferroni-adjusted Dunn)") +
			theme_classic() +
			theme(
				plot.subtitle = element_text(size = 9),
				legend.position = "none",
				aspect.ratio = 0.95
			)
		
		print(p)
		
		ggsave(
			filename = paste0("outputs/plots/final/", pc, "_by_", factor, ".png"),  # Dynamic filename
			plot = p,
			width = 6,
			height = 5,
			units = "in",
			dpi = 300
		)
	}
	
	# ── example calls
	
	# one factor across one PC axis
	print(kw_boxplot(meta_pc, pc = "PC1", factor = "Catchment"))
	
	# Each standalone factor across all PC axes
	walk(pcs, ~ print(kw_boxplot(meta_pc, pc = .x, factor = "Catchment")))
	walk(pcs, ~ print(kw_boxplot(meta_pc, pc = .x, factor = "year")))
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
}

{ # Catchment-level metal plot ----
	## pick easily recognisable symbols
	shape_vals <- c(
		"Carnon" = 16,
		"Red"    = 17,
		"Cober"  = 15,
		"Hayle"  = 18
	)
	
	# ----- 2. build the distance matrix on the four PCs
	# PCs are already scaled, so Euclidean distance is fine
	dist_mat <- dist(meta_pc %>% select(PC1:PC4), method = "euclidean")
	
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
			  x = "nMDS axis 1", y = "nMDS axis 2")
	p + ggforce::geom_mark_hull(aes(fill = Catchment),
										 concavity = 5, linetype = 0,
										 alpha = 0.08)
}


