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
		
		samp <- meta |>
			group_by(Catchment) |>
			slice_sample(n = min_n) |>
			ungroup()
		
		pca <- prcomp(samp[, metal_cols], center = TRUE, scale. = TRUE)
		
		rot_ok <- align_signs(pca$rotation[, 1:keep], ref_rot)
		
		as_tibble(rot_ok, rownames = "Metal") |>
			pivot_longer(-Metal, names_to = "PC", values_to = "loading") |>
			mutate(iter = seed)
	}
	
	# 999 balanced PCAs
	plan(multisession, workers = availableCores() - 1)
	boot_load <- future_map_dfr(1:999, function(iter) {
		if (iter %% 50 == 0) cat("Completed", iter, "of 999 iterations\n")
		one_draw(iter)
	},
	.options = furrr_options(seed = TRUE)
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

{ # validate naive bact PERMANOVA (uneven sites per river) (cross-catchment) ----
	for (i in 4) { # ready for more tax. levels if needed
		min_n <- bact_data[[i]] %>% count(Catchment) %>% pull(n) %>% min()
		
		## function that does one balanced draw + PERMANOVA
		one_run <- function(n_sub = min_n) {
			
			# 1. draw n_sub rows from every Catchment
			samp <- bact_data[[i]] %>%
				group_by(Catchment) %>%
				slice_sample(n = n_sub) %>%
				ungroup()
			
			# 3. PERMANOVA with the full formula
			mod <- adonis2(
				samp[, 21:ncol(samp)] ~
					Catchment + Season + DTM.km + Veg_Surround + Alt_Class + Temp + pH,
				data = samp,
				by = "terms",
				method = "bray",
				permutations = 999,
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
		val_overall_bact <- out %>%
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
