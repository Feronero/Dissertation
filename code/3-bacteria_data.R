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
				-Alt_class,
				-EC,
				-Site_code,
				-Site_code_2,
				-Site_time_code,
				-Poll_class
			) %>%
			relocate( Year, .after = "Location")
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

{ # Calculate median qPCR Cq Values (from the two technical replicates each) ----

	qPCR_data <- qPCR_data %>%
		group_by(Sample) %>%
		summarise(
			Cq_Median = median(Cq, na.rm = TRUE),
			Cq_IQR = IQR(Cq, na.rm = TRUE)
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
	# add "Season" and "Batch" factor columns
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- bact_data[[i]] %>%
				mutate(
					Catchment = factor(Catchment),
					Season = case_when(
						Month %in% "04" ~ "Spring",
						Month %in% "07" ~ "Summer",
						Month %in% "10" ~ "Autumn",
						Month %in% "01" ~ "Winter"
					),
					Season = factor(
						Season,
						levels = c("Spring", "Summer", "Autumn", "Winter")
					),
					Batch = factor(interaction(Site, Season, drop = TRUE)),
				) %>%
				relocate(Season, .after = "Location") %>%
				relocate(Batch, .after = "Catchment")
		}
	# convert "Catchment" and "Site to type factor
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- bact_data[[i]] %>%
				mutate(
					Catchment = factor(Catchment),
					Site = factor(Site)
				)
		}
}

{ # Rename and clean up Vegetation Labels ----
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- bact_data[[i]] %>%
				rename(Veg_Surround = Vegetation..surroundings) %>%
				mutate(Veg_Surround = factor(
					case_when(
						grepl("wooded", Veg_Surround, ignore.case = TRUE) ~ "wooded",
						grepl("open", Veg_Surround, ignore.case = TRUE) ~ "open",
						grepl("urban", Veg_Surround, ignore.case = TRUE) ~ "urban",
						grepl("overhanging", Veg_Surround, ignore.case = TRUE) ~ "wooded",
						TRUE ~ "other"
					)
				))
		}
	# View the updated dataset
		glimpse(bact_data[[1]])
}

{ # create little selection helper function ----
		select_after <- function(df, colname) {
			start <- which(names(df) == colname) + 1
			df[, start:ncol(df)]
		}
}

{ # remove un-necessary columns ----
	for (i in 1:keyvals$n_levels) {
		bact_data[[i]] <- bact_data[[i]] %>%
			select(
				Batch,
				Catchment,
				Sample,
				Site,
				Season,
				Altitude,
				DTM.km,
				Veg_Surround,
				Temp,
				pH,
				Cq_Median,
				Cq_IQR,
				richness,
				shannon,
				simpson,
				invsimpson:last_col()
			)
	}
}

{ # screen for co-linear continuous vars (only DTM and Alt are strong + sig) ----
	# only continuous explanatory vars
		vars <- c("Altitude", "DTM.km", "Temp", "pH")
	# Compute correlation coefficients and raw p-values
		res <- expand.grid(x = vars, y = vars, stringsAsFactors = FALSE) %>%
		  filter(x < y) %>%
		  mutate(
		    test = map2(x, y, ~ cor.test(
		    	bact_data[[1]][[.x]], bact_data[[1]][[.y]],
		    	method = "spearman", exact = FALSE
		    )),
		    rho = map_dbl(test, "estimate"),
		    p_raw = map_dbl(test, "p.value"),
		    p_adj = p.adjust(p_raw, method = "holm")
		  )
	# Holm-adjusted p-values
		res <- res %>%
		  mutate(p_adj = p.adjust(p_raw, method = "holm"))
	# Convert to matrices
		rho_mat <- matrix(NA, length(vars), length(vars), dimnames = list(vars, vars))
		p_raw_mat <- rho_mat; p_adj_mat <- rho_mat
		for (i in seq_len(nrow(res))) {
		  xi <- res$x[i]; yi <- res$y[i]
		  rho_mat[xi, yi] <- rho_mat[yi, xi] <- res$rho[i]
		  p_raw_mat[xi, yi] <- p_raw_mat[yi, xi] <- res$p_raw[i]
		  p_adj_mat[xi, yi] <- p_adj_mat[yi, xi] <- res$p_adj[i]
		}
		diag(rho_mat) <- 1
		diag(p_raw_mat) <- NA
		diag(p_adj_mat) <- NA
	# Display results
		list(correlations = rho_mat,
		     p_values = p_raw_mat,
		     p_adjusted = p_adj_mat)
}

cor.test(bact_data[[1]]$DTM.km, bact_data[[1]]$Altitude, method = "spearman", exact = FALSE)

{ # Merge DTM and Altitude into single shared PC axis ----
	# PCA combining DTM.km and Altitude
		pca_result <- prcomp(
			bact_data[[1]] %>%
				select(DTM.km, Altitude),
			scale. = TRUE,
			center = TRUE
		)
		summary(pca_result)$importance %>%
			print()
		pca_result$rotation %>%
			print()
	# Replace DTM.km and Altitude with PC1 ("River_Gradient")
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] <- bact_data[[i]] %>%
				mutate(River_Gradient = -pca_result$x[, 1]) %>%
				select(-DTM.km, -Altitude) %>%
				relocate(River_Gradient, .after = "Season")
		}
}

{ # re-screen for co-linear continuous vars after combination (no further corrs) ----
	# only continuous explanatory vars
		vars <- c("River_Gradient", "Temp", "pH")
	# Compute correlation coefficients and raw p-values
		res <- expand.grid(x = vars, y = vars, stringsAsFactors = FALSE) %>%
		  filter(x < y) %>%
		  mutate(
		    test = map2(x, y, ~ cor.test(
		    	bact_data[[1]][[.x]], bact_data[[1]][[.y]],
		    	method = "spearman", exact = FALSE
		    )),
		    rho = map_dbl(test, "estimate"),
		    p_raw = map_dbl(test, "p.value"),
		    p_adj = p.adjust(p_raw, method = "holm")
		  )
	# Holm-adjusted p-values
		res <- res %>%
		  mutate(p_adj = p.adjust(p_raw, method = "holm"))
	# Convert to matrices
		rho_mat <- matrix(NA, length(vars), length(vars), dimnames = list(vars, vars))
		p_raw_mat <- rho_mat; p_adj_mat <- rho_mat
		for (i in seq_len(nrow(res))) {
		  xi <- res$x[i]; yi <- res$y[i]
		  rho_mat[xi, yi] <- rho_mat[yi, xi] <- res$rho[i]
		  p_raw_mat[xi, yi] <- p_raw_mat[yi, xi] <- res$p_raw[i]
		  p_adj_mat[xi, yi] <- p_adj_mat[yi, xi] <- res$p_adj[i]
		}
		diag(rho_mat) <- 1
		diag(p_raw_mat) <- NA
		diag(p_adj_mat) <- NA
	# Display results
		list(correlations = rho_mat,
		     p_values = p_raw_mat,
		     p_adjusted = p_adj_mat)
}

{ # Test VIF (multivariate co-linearity) (Temp moderately inflated VIF, maybe drop) ----
	lm(
		River_Gradient ~ 
			Catchment + Season + Veg_Surround + pH + Temp,
		data = bact_data[[1]]
	) %>%
	check_collinearity()
}

{ # Test dependent variables for co-linearity ----
	test <- list()
		for (i in 1:keyvals$n_levels) {
			bact_data[[i]] %>%
				select(invsimpson, Cq_Median, richness) %>%
				drop_na() %>%
				cor(method = "spearman") %>%
				print()
			pca_res <- bact_data[[i]] %>%
				drop_na(Cq_Median) %>%
				select(invsimpson, richness) %>%
				scale() %>%
				prcomp()
			pca_res %>%
				summary() %>%
				print()
			test[[i]] <- bact_data[[i]]
			test[[i]] <- test[[i]] %>%
				drop_na(Cq_Median) %>%
				mutate(
					PC1 = pca_res$x[, 1],
					PC2 = pca_res$x[, 2],
					pH = as.numeric(scale(pH)),
					Temp = as.numeric(scale(Temp)),
					River_Gradient = as.numeric(scale(River_Gradient))
				) %>%
				relocate(
					PC1, PC2,
					.after = pH
				)
		}
}
	
{ # export and clean-up temp objects ----

	# export the level 4 bact_data as a .csv
		write_csv(
			x = bact_data[[1]],
			file = "outputs/tables/bact_data.csv"
		)
	# remove individual dataframes
		rm(
			a_div, bact_meta, PCR_data, qPCR_data, p_adj_mat, p_raw_mat, pca_res, pca_result, res, rho_mat
		)
}

{ # Check subsampling depth ----
	level = 1
	rowSums(bact_data[[level]] %>% select_after("invsimpson"))
}
