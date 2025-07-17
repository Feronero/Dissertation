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

{ # collapse the data to one level per site/season ----
	bact_collapsed <- list()
	for (i in 1:keyvals$n_levels) {
		bact_collapsed[[i]] <- as.data.frame(bact_data[[i]] %>%
			group_by(Batch) %>%
			summarise(
				Catchment = first(Catchment),
				Site = first(Site),
				Season = first(Season),
				Altitude = first(Altitude),
				DTM.km = first(DTM.km),
				Veg_Surround = first(Veg_Surround),
				Temp = first(Temp),
				pH = first(pH),
				Cq_Median = median(Cq_Median, na.rm = TRUE),
				Cq_IQR = median(Cq_IQR, na.rm = TRUE),
				richness = median(richness, na.rm = TRUE),
				shannon = median(shannon, na.rm = TRUE),
				simpson = median(simpson, na.rm = TRUE),
				invsimpson = median(invsimpson, na.rm = TRUE),
				across(
					names(bact_data[[i]] %>% select_after("invsimpson")),
					\(x) median(x, na.rm = TRUE)
				)
			) %>%
			ungroup()
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
		    	bact_collapsed[[1]][[.x]], bact_collapsed[[1]][[.y]],
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

{ # Merge DTM and Altitude into single shared PC axis ----
	# PCA combining DTM.km and Altitude
		pca_result <- prcomp(
			bact_collapsed[[1]] %>%
				select(DTM.km, Altitude),
			scale. = TRUE,
			center = TRUE
		)
		summary(pca_result)$importance
		pca_result$rotation
	# Replace DTM.km and Altitude with PC1 ("River_Gradient")
		for (i in 1:keyvals$n_levels) {
			bact_collapsed[[i]] <- bact_collapsed[[i]] %>%
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
		    	bact_collapsed[[1]][[.x]], bact_collapsed[[1]][[.y]],
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

{ # Test VIF (multivariate co-linearity) (Temp moderately infalted VIF, maybe drop) ----
	vif_out <- vif(lm(
		River_Gradient ~ 
			Catchment + Season + Veg_Surround + pH + Temp,
		data = bact_collapsed[[1]]
	))
	print(vif_out)
}

{ # Compare PERMANOVA performance with + without Temp ----
	# Fit both PERMANOVA models
		mod_with <- adonis2(
			bact_collapsed[[4]] %>% select_after("invsimpson") ~
				Catchment + Season + River_Gradient + Veg_Surround + pH + Temp,
			data = bact_collapsed[[4]],
			method = "bray",
			by = "margin",
			strata = bact_collapsed[[4]]$Site
		)
		mod_no <- adonis2(
			bact_collapsed[[4]] %>% select_after("invsimpson") ~
				Catchment + Season + River_Gradient + Veg_Surround + pH,
			data = bact_collapsed[[4]],
			method = "bray",
			by = "margin",
			strata = bact_collapsed[[4]]$Site
		)
	# Compute AICc for each model
		aicc_with <- AICc_permanova2(mod_with)$AICc
		aicc_no   <- AICc_permanova2(mod_no)$AICc
	# Extract total R2 and adjusted R2
		getR2 <- function(mod) {
		  ss      <- mod$SumOfSqs
		  totalSS <- sum(ss)
		  R2_tot  <- (sum(ss[-length(ss)]) / totalSS)
		  n <- sum(mod$Df)
		  p <- length(ss) - 1
		  R2_adj <- 1 - (1 - R2_tot) * ((n - 1) / (n - p - 1))
		  return(c(R2 = R2_tot, Adj_R2 = R2_adj))
		}
		r2_with <- getR2(mod_with)
		r2_no   <- getR2(mod_no)
	# Combine results into a comparison table
		out <- tibble(
		  Model     = c("With Temp", "Without Temp"),
		  R2        = c(r2_with["R2"], r2_no["R2"]),
		  Adj_R2    = c(r2_with["Adj_R2"], r2_no["Adj_R2"]),
		  AICc      = c(aicc_with, aicc_no)
		)
		print(out)
	# clean up
		rm(mod_with)
}

{ # View Final PERMANOVA results ---
	mod_no
}

{ # Catchment-level nMDS plot ----
	# pick easily recognisable symbols
	shape_vals <- c(
		"Spring" = 16,
		"Summer"    = 17,
		"Autumn"  = 15,
		"Winter"  = 18
	)
	
	# run non-metric MDS (two dimensions)
	set.seed(123)
	fit <- metaMDS(
		comm = bact_collapsed[[4]] %>% select_after("invsimpson"),
		distance = "bray",
		k = 2,
		trymax = 100
	)
	
	# combine ordination scores with the metadata
	plot_df <- bact_collapsed[[4]] %>%
		mutate(
			MDS1 = fit$points[, 1],
			MDS2 = fit$points[, 2]
		)
	
	# 5. draw the map
	p <- ggplot(
		plot_df,
			aes(
				MDS1, MDS2,
				colour = Season,
				shape = Season
			)
		) +
		geom_point(size = 3, alpha = 0.8) +
		scale_shape_manual(values = shape_vals) +
		theme_classic() +
		labs(
			title = NULL,
			x = "nMDS Axis 1",
			y = "nMDS Axis 2"
		)
	p <- p +
		geom_mark_hull(
			aes(fill = Season),
			concavity = 5,
			linetype = 0,
			alpha = 0.08
		)
	
	ggsave(
		filename = paste0("outputs/plots/final/bacteria/Seasons.png"),
		plot = p,
		width = 9,
		height = 5,
		units = "in",
		dpi = 1000
	)
	p
}

{ # db-RDA of each factor individually v1 ----
	# Pull community matrix
		comm <- bact_collapsed[[4]] %>% select_after("invsimpson")
	# Store db-RDA plots
		db_rda_plots <- list()
	# Factors to test
		factors <- c("Catchment", "Season", "River_Gradient", "Veg_Surround")
	# Loop through each factor
		for (var in factors) {
	# If numeric, bin to factor for plotting
		  var_data <- bact_collapsed[[4]][[var]]
		  if (is.numeric(var_data)) {
		    bact_collapsed[[4]][[paste0(var, "_binned")]] <- cut(var_data, breaks = 4)
		    var_plot <- paste0(var, "_binned")
		  } else {
		    var_plot <- var
		  }
	# Model
		  form <- as.formula(paste0("comm ~ ", var))
		  mod <- capscale(form, data = bact_collapsed[[4]], distance = "bray")
	# Scores
		  scores_df <- bact_collapsed[[4]] %>%
		    mutate(
		      CAP1 = scores(mod, display = "sites")[, 1],
		      CAP2 = scores(mod, display = "sites")[, 2],
		      group = .data[[var_plot]]
		    )
	# Hulls
		  hulls_df <- scores_df %>%
		    group_by(group) %>%
		    slice(chull(CAP1, CAP2))
		  
	# Plot
		  p <- ggplot(scores_df, aes(CAP1, CAP2)) +
		    geom_polygon(
		      data = hulls_df,
		      aes(fill = group, group = group),
		      alpha = 0.08,
		      colour = NA
		    ) +
		    geom_point(
		      aes(colour = group, shape = group),
		      size = 3,
		      alpha = 0.8
		    ) +
		    theme_classic() +
		    labs(
		      x = "db-RDA Axis 1",
		      y = "db-RDA Axis 2",
		      title = paste("db-RDA:", var)
		    ) +
		    guides(fill = "none")  # hide fill legend if it duplicates colour
		  
	# Save
		  ggsave(
		    filename = paste0("outputs/plots/final/bacteria/dbRDA_", var, ".png"),
		    plot = p,
		    width = 9,
		    height = 5,
		    units = "in",
		    dpi = 1000
		  )
		  
		  db_rda_plots[[var]] <- p
		}
}

{ # View v1 db-RDA plots ----
		db_rda_plots$Catchment
		db_rda_plots$Season
		db_rda_plots$Veg_Surround
		db_rda_plots$River_Gradient
}

{ # db-RDA of each factor individually v2 ----
		alpha_max <- 0.4
		alpha_min <- 0.2
		density_cutoff_quantile <- 0.6
		
		db_rda_plots <- list()
		factors <- c("Catchment", "Season", "River_Gradient", "Veg_Surround")
		
		for (var in factors) {
		  var_data <- bact_collapsed[[4]][[var]]
		  group_var <- if (is.numeric(var_data))
		    cut(var_data, breaks = 4)
		  else var_data
		
		  mod <- capscale(as.formula(paste0("comm ~ ", var)),
		                 data = bact_collapsed[[4]], distance = "bray")
		
		  scores_df <- bact_collapsed[[4]] %>% transmute(
		    CAP1 = scores(mod, display = "sites")[,1],
		    CAP2 = scores(mod, display = "sites")[,2],
		    group = group_var
		  )
		
		  group_levels <- sort(unique(as.character(scores_df$group)))
		  n_groups <- length(group_levels)
		
		  shape_vals <- setNames(21:(20+n_groups), group_levels)
		  group_colors <- setNames(
		    hcl(seq(15, 375, length.out = n_groups + 1)[-1], 100, 65),
		    group_levels
		  )
		
		  blur_layers <- scores_df %>%
		    group_split(group, .keep = TRUE) %>%
		    map(function(df) {
		      this_group <- unique(as.character(df$group))
		      ddata <- ggplot_build(
		        ggplot(df, aes(CAP1, CAP2)) +
		        stat_density_2d(geom = "tile", contour = FALSE, n = 150)
		      )$data[[1]]
		
		      max_d <- max(ddata$density, na.rm = TRUE)
		      cutoff <- quantile(ddata$density, density_cutoff_quantile, na.rm = TRUE)
		
		      stat_density_2d(
		        data = df,
		        aes(
		          x = CAP1, y = CAP2,
		          alpha = ifelse(
		            after_stat(density) < cutoff, 0,
		            alpha_min + (after_stat(density)-cutoff)/(max_d-cutoff)*(alpha_max-alpha_min)
		          ),
		          fill = this_group
		        ),
		        geom = "tile", contour = FALSE, n = 600,
		        inherit.aes = FALSE, show.legend = FALSE
		      )
		    })
		
		  p <- ggplot(scores_df, aes(CAP1, CAP2)) +
		    blur_layers +
		    geom_point(aes(fill = group, shape = group),
		               colour = "black", stroke = 0,
		               size = 3, alpha = 1) +
		    scale_shape_manual(values = shape_vals, name = NULL) +
		    scale_fill_manual(values = group_colors, name = NULL) +
		    scale_alpha_identity() +
		    guides(
		      fill = guide_legend(override.aes = list(
		        shape = unname(shape_vals),
		        fill = unname(group_colors),
		        colour = "black",
		        stroke = 0,
		        size = 3,
		        alpha = 1
		      )),
		      shape = "none"
		    ) +
		    theme_classic() +
		    labs(x = "db-RDA Axis 1", y = "db-RDA Axis 2")
		
		  ggsave(
		    paste0("outputs/plots/final/bacteria/dbRDA_", var, ".png"),
		    plot = p, width = 9, height = 5, units = "in", dpi = 1000
		  )
		
		  db_rda_plots[[var]] <- p
		}
}

{ # View updated db-RDA plots ----
	db_rda_plots$Season
	db_rda_plots$Catchment
	db_rda_plots$River_Gradient
	db_rda_plots$Veg_Surround
}
	
{ # export and clean-up temp objects ----

	# export the level 4 bact_data as a .csv
		write_csv(
			x = bact_collapsed[[1]],
			file = "outputs/tables/bact_collapsed.csv"
		)
	# remove individual dataframes
		rm(
			a_div, bact_meta, PCR_data, qPCR_data
		)
	
}

{ # Check subsampling depth ----
	level = 1
	rowSums(bact_data[[level]] %>% select_after("invsimpson"))
}
