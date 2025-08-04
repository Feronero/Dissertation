# This script allows all the long validation PERMANOVAs and PCAs to be run in the background
# Requires packages to be loaded in explicitly

{ # Packages ----
	library(furrr)
	library(future)
#	library(agricolae)
	library(vegan)
	library(tidyverse)
	library(permute)
#	library(FactoMineR)
#	library(factoextra)
#	library(writexl)
	#	library(MASS)			# conflicts with dplyr
#	library(viridisLite)
	#	library(fitdistrplus)	# requires MASS
#	library(mgcv)
#	library(car)      # For Levene's test
#	library(rstatix)  # For pairwise comparisons
#	library(lubridate)
	library(ggforce)
	library(ggpubr)
}

{ # Validate unbalanced PCA loadings ----
	# calculate mininumum n metal samples per catchment
		min_n <- meta %>%
			count(Catchment) %>%
			pull(n) %>%
			min()
	# reference rotation – enforce metal_cols order
	ref_rot <- naive_loadings_df %>%
		arrange(match(Metal, metal_cols)) %>%   # <<< new line
		column_to_rownames("Metal") %>%
		as.matrix()
	
	# helper unchanged
align_signs <- function(new_rot, ref_rot) {
  ref <- ref_rot[rownames(new_rot), , drop = FALSE]
  cor_mat <- cor(new_rot, ref)
  new_order <- apply(cor_mat, 2, which.max)
  
  # Ensure no duplicates in new_order
  if (any(duplicated(new_order))) {
    warning("Duplicate column matches detected. Using greedy matching.")
    # Alternative: Assign unmatched PCs to remaining columns
    unused <- setdiff(seq_len(ncol(new_rot)), new_order)
    for (i in which(duplicated(new_order))) {
      new_order[i] <- unused[1]
      unused <- unused[-1]
    }
  }
  
  new_rot <- new_rot[, new_order, drop = FALSE]
  colnames(new_rot) <- colnames(ref)  # Preserve reference names
  cor_sign <- sign(diag(cor(new_rot, ref)))
  sweep(new_rot, 2, cor_sign, `*`)
}
	
	# one balanced bootstrap draw
	one_draw <- function(seed) {
		
		samp <- meta %>%
			group_by(Catchment) %>%
			slice_sample(n = min_n) |>
			ungroup()
		
		pca <- prcomp(samp[, metal_cols], center = TRUE, scale. = TRUE)
		
		rot_ok <- align_signs(pca$rotation[, 1:keep], ref_rot)
		
		as_tibble(rot_ok, rownames = "Metal") |>
			pivot_longer(-Metal, names_to = "PC", values_to = "loading") |>
			mutate(iter = seed)
	}
	
	# 999 balanced PCAs
	plan(multisession, workers = 8)
	boot_load <- future_map_dfr(
		1:999, function(iter) {
		one_draw(iter)
	},
	.options = furrr_options(seed = TRUE),
	.progress = TRUE
	)
	
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
	
	val_naive_PCA <- compare %>%
		group_by(PC) %>%
		summarise(correlation = cor(median_loading, orig_loading),
					 rmse        = sqrt(mean((orig_loading - median_loading)^2)),
					 pct_within  = mean(within_99) * 100,
					 .groups = "drop")
	
#	write_csv(compare, "outputs/tables/loading_comparison_detail.csv")
#	write_csv(val_naive_PCA,  "outputs/tables/val_naive_PCA.csv")
}

{ # Entirely bootstrapped PCA ----
	
	# Set up parallel processing
	plan(multisession, workers = 8)
	
	# Function to perform PCA on a balanced subsample and return scores and loadings
	balanced_pca_with_metrics <- function(data, min_n, pca_vars) {
	  # Subsample data to balanced design
	  subsampled_data <- data %>%
	    group_by(Catchment) %>%
	    sample_n(size = min_n, replace = FALSE) %>%
	    ungroup()
	  
	  # Store the sample IDs that were selected
	  selected_samples <- subsampled_data %>% select(notation, Catchment, Site, year, Season, month)
	  
	  # Perform PCA on metal concentrations
	  pca_result <- subsampled_data %>%
	    select(all_of(pca_vars)) %>%
	    prcomp(center = TRUE, scale. = TRUE)
	  
	  # Get PC scores for the first 3 components
	  pc_scores <- as.data.frame(pca_result$x[, 1:3])
	  names(pc_scores) <- c("PC1", "PC2", "PC3")
	  
	  # Get loadings for the first 3 components
	  pc_loadings <- as.data.frame(pca_result$rotation[, 1:3])
	  names(pc_loadings) <- c("PC1_loading", "PC2_loading", "PC3_loading")
	  pc_loadings$metal <- rownames(pc_loadings)
	  
	  # Get variance explained
	  var_explained <- summary(pca_result)$importance[2, 1:3] # Proportion of variance
	  var_df <- data.frame(
	    PC = c("PC1", "PC2", "PC3"),
	    var_explained = var_explained
	  )
	  
	  # Combine results
	  list(
	    samples = cbind(selected_samples, pc_scores, iteration = attr(data, "iteration")),
	    loadings = cbind(pc_loadings, iteration = attr(data, "iteration")),
	    variance = cbind(var_df, iteration = attr(data, "iteration"))
	  )
	}
	
	# Main analysis function to get aggregated metrics
	run_balanced_pca_with_metrics <- function(metal_wide, n_iterations = 8000) {
	  # Define which columns are metal concentrations
	  pca_vars <- c("Cadmium", "Iron", "Nickel", "Copper", 
	                "Manganese", "Zinc", "Magnesium")
	  
	  # Calculate minimum sample size per catchment
	  min_n <- metal_wide %>%
	    count(Catchment) %>%
	    pull(n) %>%
	    min()
	  
	  # Create iteration counter attribute
	  metal_wide_list <- map(1:n_iterations, ~{
	    df <- metal_wide
	    attr(df, "iteration") <- .x
	    df
	  })
	  
	  # Run parallel balanced PCAs and collect results
	  future_map(
	    metal_wide_list,
	    ~balanced_pca_with_metrics(.x, min_n, pca_vars),
	    .options = furrr_options(seed = TRUE),
	    .progress = TRUE
	  )
	}
	
	# Run the analysis to get all PCA metrics
	all_pca_results <- run_balanced_pca_with_metrics(metal_wide)
	
	# Extract and combine all components
	all_pc_scores <- map_dfr(all_pca_results, ~.x$samples)
	all_loadings <- map_dfr(all_pca_results, ~.x$loadings)
	all_variance <- map_dfr(all_pca_results, ~.x$variance)
	
	# Calculate aggregated PC scores for each sample
	aggregated_scores <- all_pc_scores %>%
	  group_by(notation, Catchment, Site, year, Season, month) %>%
	  summarise(
	    PC1_median = median(PC1, na.rm = TRUE),
	    PC1_iqr = IQR(PC1, na.rm = TRUE),
	    PC2_median = median(PC2, na.rm = TRUE),
	    PC2_iqr = IQR(PC2, na.rm = TRUE),
	    PC3_median = median(PC3, na.rm = TRUE),
	    PC3_iqr = IQR(PC3, na.rm = TRUE),
	    n_included = n(),  # How many times each sample was included
	    .groups = "drop"
	  )
	
	# Calculate median loadings and IQR for each metal in each PC
	aggregated_loadings <- all_loadings %>%
	  group_by(metal) %>%
	  summarise(
	    PC1_loading_median = median(PC1_loading, na.rm = TRUE),
	    PC1_loading_iqr = IQR(PC1_loading, na.rm = TRUE),
	    PC2_loading_median = median(PC2_loading, na.rm = TRUE),
	    PC2_loading_iqr = IQR(PC2_loading, na.rm = TRUE),
	    PC3_loading_median = median(PC3_loading, na.rm = TRUE),
	    PC3_loading_iqr = IQR(PC3_loading, na.rm = TRUE),
	    n_included = n(),  # How many times each metal was included
	    .groups = "drop"
	  )
	
	# Calculate median variance explained and IQR for each PC
	aggregated_variance <- all_variance %>%
	  group_by(PC) %>%
	  summarise(
	    var_explained_median = median(var_explained, na.rm = TRUE),
	    var_explained_iqr = IQR(var_explained, na.rm = TRUE),
	    n_included = n(),  # How many times each PC was included
	    .groups = "drop"
	  )
	
	# Join the aggregated scores back to the original data
	metal_wide_with_scores <- metal_wide %>%
	  left_join(aggregated_scores, by = c("notation", "Catchment", "Site", "year", "Season", "month"))
	
	# View the results
	bootstrap_PCA_res <- list(
	  samples = metal_wide_with_scores,
	  loadings = aggregated_loadings,
	  variance = aggregated_variance
	)
	
}
saveRDS(bootstrap_PCA_res, "outputs/backups/bootstrap_PCA_res.rds")
1 <- readRDS("outputs/backups/bootstrap_PCA_res.rds")

# add the results to metal_data_complete
	metal_data_complete <- bootstrap_PCA_res$samples %>%
		select(c(1:13, 14, 16, 18))
saveRDS(metal_data_complete, "outputs/backups/metal_data_complete.rds")
metal_data_complete <- readRDS("outputs/backups/metal_data_complete.rds")

{ # Sub-Sampled cross-catch metal PERMANOVAs ----
	# even samps per catchment
	min_n <- metal_data_complete %>%
	    count(Catchment) %>%
	    pull(n) %>%
	    min()	
	# function for one balanced iteration
		one_run <- function() {
	# randomly subsample from each catchment
			samp <- metal_data_complete %>%
				group_by(Catchment) %>%
				slice_sample(n = min_n) %>%
				ungroup()
	# establish permutation restrictions
			perm <- how(
				blocks = samp$Catchment,
				nperm = 999
			)
	# run one PERMANOVA iteration
			mod <- adonis2(
				samp %>% select(PC1_median, PC2_median, PC3_median) ~
					Catchment + year + Season + Site,
				data = samp,
				permutations = perm,
				method = "euclidean",
				by = "margin",
				parallel = 2
			)
	# save key output
			tibble(
				term = rownames(mod),
				R2   = mod$R2,
				p    = mod$`Pr(>F)`
			)
		}
	# set initial seed for reproducibility
		set.seed(42)
	# run all iterations via multi-core
		plan(multisession, workers = 4)
		out <- future_map_dfr(
			1:999,
			~ one_run(),
			.options = furrr_options(seed = TRUE),
			.progress = TRUE
		)
	# save summary stats of all runs
		metal_PERM <- out %>%
			group_by(term) %>%
			summarise(
				median_R2 = median(R2),
				lo_R2 = quantile(R2, .025),
				hi_R2 = quantile(R2, .975),
				`%p<0.05` = mean(p < 0.05) * 100,
				`%p<0.01` = mean(p < 0.01) * 100,
				`%p<0.001` = mean(p < 0.001) * 100,
				.groups = "drop"
			)
}
metal_PERM # balanced sites per catch, VERY NOT balanced sampls per site
saveRDS(metal_PERM, "outputs/backups/metal_PERM.rds")
metal_PERM <- readRDS("outputs/backups/metal_PERM.rds")

# massive imbalance in samples per site
min_per_site <- metal_data_complete %>%
			count(Catchment, Site) %>%
			group_by(Catchment) %>%
			summarise(m = min(n)) %>%
			pull(m) %>% min()

{ # validate unbalanced PC-score (metal) PERMANOVA using naive PC scores  ----
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
		dist_mat <- dist(samp %>% select(all_of(paste0("PC", 1:keep))))
		
		# 3. PERMANOVA with the full formula
		mod <- adonis2(
			dist_mat ~ Catchment + Site %in% Catchment + year %in% Catchment + Season,
			data = samp,
			by = "terms",
			permutations = 999,
			parallel = 4
		)
		
		# 4. keep just the rows you care about
		tibble(term = rownames(mod)[1:keep],
				 R2   = mod$R2[1:keep],
				 p    = mod$`Pr(>F)`[1:keep])
	}
	
	## 999 balanced runs
	set.seed(42)
	plan(multisession, workers = 3)
	out <- future_map_dfr(1:999, function(iter) {
		if (iter %% 50 == 0) cat("Completed", iter, "of 999 iterations\n")
		one_run()
	})
	
	## summarise
	val_PC_PERMANOVA <- out %>%
		group_by(term) %>%
		summarise(median_R2 = median(R2),
					 lo_R2     = quantile(R2, 0.025),
					 hi_R2     = quantile(R2, 0.975),
					 pct_p_lt_0.05 = mean(p < 0.05) * 100,
					 .groups = "drop")
}

{ # validate naive bact-metal PERMANOVA (each catchment) ----
		val_metalbact_PERM_percatch <- list()
	
		for (i in 4) { # ready for other tax. levels if needed
			for (j in keyvals$catchments) {
			
	# create catchment-specific subset + count min samples per site
				catch_subset <- master[[i]] %>% filter(Catchment == j)
				min_n <- catch_subset %>% count(Site) %>% pull(n) %>% min()
			
	# function that does one balanced draw + PERMANOVA
				one_run <- function(n_sub = min_n) {
				
	# draw n_sub rows from every Catchment
					samp <- catch_subset %>%
						group_by(Site) %>%
						slice_sample(n = n_sub) %>%
						ungroup()
				
	# run subsampled PERMANOVA
					mod <- adonis2(
						samp[, 24:ncol(samp)] ~
							Site + Season + DTM.km + Veg_Surround + Alt_Class + Temp + pH + PC1 + PC2 + PC3,
						data = samp,
						by = "terms",
						method = "bray",
						permutations = 999,
						parallel = 4
					)
	# results of single run
					tibble(
						term = rownames(mod),
						R2   = mod$R2,
						p    = mod$`Pr(>F)`
					)
				}
	# start 999 balanced runs
				set.seed(42)
				plan(multisession, workers = 3)
				out <- map_dfr(1:999, function(iter) {
					if (iter %% 50 == 0) cat("Completed", iter, "of 999 iterations\n")
					one_run()
				})
			
	# summarise overall results
				val_metalbact_PERM_percatch[[j]] <- out %>%
					group_by(term) %>%
					summarise(
						median_R2 = median(R2),
						# these give the 95% interval
						lo_R2     = quantile(R2, 0.025),
						hi_R2     = quantile(R2, 0.975),
						pct_p_lt_0.05 = mean(p < 0.05) * 100,
						.groups = "drop"
					)
		}
	}
}

{ # validate naive bact-metal PERMANOVA (cross-catchment) ----
	# uses CarnonRed_collapsed - minus Hayle and Cober due to insufficient sample size to accept the full formula for verification
	for (i in 4) { # ready for more tax. levels if needed
		min_n <- CarnonRed_collapsed[[i]] %>% dplyr::count(Catchment) %>% pull(n) %>% min()
		
		## function that does one balanced draw + PERMANOVA
		one_run <- function(n_sub = min_n) {
			
			# draw n_sub rows from every Catchment
			samp <- CarnonRed_collapsed[[i]] %>%
				group_by(Catchment) %>%
				slice_sample(n = n_sub) %>%
				ungroup()
			
			samp_OTUs <- samp %>% select_after("PC3")
			
			# 3. PERMANOVA with the full formula
			mod <- adonis2(
				samp_OTUs ~
					Catchment + Season + PC2,
				data = samp,
				by = "terms",
				method = "bray",
				permutations = 999,
				strata = samp$Site,
				parallel = 4
			)
			
			# 4. keep just the rows you care about
			tibble(
				term = rownames(mod),
				R2   = mod$R2,
				p    = mod$`Pr(>F)`
			)
		}
		
		## 999 balanced runs
		set.seed(42)
		plan(multisession, workers = 3)
		out <- future_map_dfr(1:999, function(iter) {
			if (iter %% 50 == 0) cat("Completed", iter, "of 999 iterations\n")
			one_run()
		})
		
		## summarise
		val_crsscatch_metal_PC <- out %>%
			group_by(term) %>%
			summarise(median_R2 = median(R2),
						 # these give the 95% interval
						 lo_R2     = quantile(R2, 0.025),
						 hi_R2     = quantile(R2, 0.975),
						 pct_p_lt_0.05 = mean(p < 0.05) * 100,
						 .groups = "drop")
	}
}

{ # validate naive bact-metal PERMANOVA (each catchment) ----
		val_metalbact_PERM_percatch <- list()
	
		for (i in 4) { # ready for other tax. levels if needed
			for (j in keyvals$catchments) {
			
	# create catchment-specific subset + count min samples per site
				catch_subset <- master[[i]] %>% filter(Catchment == j)
				min_n <- catch_subset %>% count(Site) %>% pull(n) %>% min()
			
	# function that does one balanced draw + PERMANOVA
				one_run <- function(n_sub = min_n) {
				
	# draw n_sub rows from every Catchment
					samp <- catch_subset %>%
						group_by(Site) %>%
						slice_sample(n = n_sub) %>%
						ungroup()
				
	# run subsampled PERMANOVA
					mod <- adonis2(
						samp[, 24:ncol(samp)] ~
							Site + Season + DTM.km + Veg_Surround + Alt_Class + Temp + pH + PC1 + PC2 + PC3,
						data = samp,
						by = "terms",
						method = "bray",
						permutations = 999,
						parallel = 4
					)
	# results of single run
					tibble(
						term = rownames(mod),
						R2   = mod$R2,
						p    = mod$`Pr(>F)`
					)
				}
	# start 999 balanced runs
				set.seed(42)
				plan(multisession, workers = 3)
				out <- map_dfr(1:999, function(iter) {
					if (iter %% 50 == 0) cat("Completed", iter, "of 999 iterations\n")
					one_run()
				})
			
	# summarise overall results
				val_metalbact_PERM_percatch[[j]] <- out %>%
					group_by(term) %>%
					summarise(
						median_R2 = median(R2),
						# these give the 95% interval
						lo_R2     = quantile(R2, 0.025),
						hi_R2     = quantile(R2, 0.975),
						pct_p_lt_0.05 = mean(p < 0.05) * 100,
						.groups = "drop"
					)
		}
	}
}
