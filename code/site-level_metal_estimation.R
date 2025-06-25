# this pipeline:
	# directly assigns metal levels to sites with metal records in the same or adjacent months
	# then estimates the emtal levels of all other points baseed on multiple stat. tests

{ # Directly pair metal and bacterial data taken in the same or adjacent months ----
	# load & tag
	bact <- bact_data[[4]] %>% 
		mutate(across(c(Year, Month), as.integer),
				 b_date = make_date(Year, Month, 1),
				 b_id   = row_number())          # unique id per bacterial row
	
	metal <- metal_wide %>% 
		rename(Year = year, Month = month) %>% 
		mutate(across(c(Year, Month), as.integer),
				 m_date = make_date(Year, Month, 1))
	
	metal_cols <- c("Cadmium","Iron","Nickel","Copper",
						 "Manganese","Zinc","Magnesium")
	
	meta_cols <- c("SampleID","Sample","Site","Catchment","Location",
						"Year","Month","Date")
	
	species_cols <- setdiff(names(bact), c(meta_cols, "b_date", "b_id"))
	
	# collect candidate bacteria:metal matches (±1 month)
	candidates <- bact %>% 
		crossing(offset = -1:1) %>%                       # replicate rows
		mutate(match_date = b_date %m+% months(offset)) %>% 
		left_join(metal %>% select(Site, m_date, all_of(metal_cols)),
					 by = c("Site",
					 		 "match_date" = "m_date"),
					 relationship = "many-to-many") %>%      # warning gone
		filter(!if_all(all_of(metal_cols), is.na))         # keep rows with metal data
	
	# choose one value for each metal element per bacterial sample
	pick_val <- function(val, off) {
		same <- val[off == 0 & !is.na(val)]
		if(length(same))                return(same[1])
		adj  <- val[off != 0 & !is.na(val)]
		if(length(adj) == 1)            return(adj)
		if(length(adj) >  1)            return(median(adj))
		NA_real_
	}
	# list the final results and label them according to the temporal source of their metal data
	paired_df <- candidates %>% 
		group_by(b_id) %>% 
		summarise(across(all_of(meta_cols),   first),
					 across(all_of(species_cols),first),
					 across(all_of(metal_cols),  ~ pick_val(.x, offset)),
					 method = case_when(
					 	any(offset == 0)          ~ "same month",
					 	n_distinct(offset) == 1   ~ "adjacent month",
					 	TRUE                      ~ "median values"
					 ),
					 .groups = "drop") %>% 
		relocate(all_of(metal_cols), .after = invsimpson) %>%
		relocate(method, .after = all_of(metal_cols))
	
	# export the results
	write_csv(
		paired_df,
		"outputs/tables/paired_df.csv"
	)
}

{ # Per-catchment graphs of all sampling sites and all metals, separated by season ----
	# separate non-data columns
	id_cols <- c("notation", "Catchment", "Site", "year", "Season", "Month")
	# pivot data
	metal_long <- metal_wide %>%                            
		pivot_longer(
			cols      = -all_of(id_cols),
			names_to  = "Metal",
			values_to = "value"
		)
	# order season factor levels
	metal_long$Season <- factor(
		metal_long$Season,
		levels = c("Spring", "Summer", "Autumn", "Winter"))
	# helper function
	plot_catchment <- function(df) {
		# add a per-group count to tell boxes from singletons
		df <- df %>% add_count(Metal, Site, Season, name = "n")
		ggplot() +
			# draw boxes only where we have ≥3 observations
			geom_boxplot(
				data  = df %>% filter(n >= 3),
				aes(Season, value),
				fill  = "grey90",
				outlier.size = 0.5
			) +
			# draw X's for 1- or 2-point groups
			geom_point(
				data  = df %>% filter(n < 3),
				aes(Season, value),
				shape = 1,                 # an “X”
				size  = 2
			) +
			facet_grid(Metal ~ Site, scales = "free_y") +
			labs(
				title = unique(df$Catchment),
				x     = "Season",
				y     = "Concentration"
			) +
			theme_classic(base_size = 9) +
			theme(
				axis.text.x     = element_text(angle = 45, hjust = 1),
				panel.spacing.y = unit(0.4, "lines")
			)
	}
	# build the four plots
	plots <- metal_long %>%
		group_split(Catchment) %>%           # one tibble per catchment
		map(plot_catchment)                  # one ggplot per tibble
	# view
	plots[[1]]
	plots[[2]]
	plots[[3]]
	plots[[4]]
}

{ # Are there significant seasonal differences in metal levels within sites? ----
	# Reshape data to long format and filter
	metal_long <- metal_wide %>%
		pivot_longer(
			cols = c(Cadmium, Iron, Nickel, Copper, Manganese, Zinc, Magnesium),
			names_to = "Metal",
			values_to = "Concentration"
		) %>%
		mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn"))) %>%
		group_by(Site, Metal) %>%
		mutate(group_count = n()) %>%
		filter(group_count >= 3) %>%
		select(-group_count) %>%
		ungroup()
	
	# Function to test seasonal differences
	test_seasonal_differences <- function(df) {
		# Check if enough seasons exist
		seasons_with_data <- length(unique(df$Season))
		if (seasons_with_data < 2) {
			return(data.frame(
				Test = "Insufficient data",
				P.value = NA,
				Significant.Pairs = "Less than 2 seasons with data",
				stringsAsFactors = FALSE
			))
		}
		
		# Check if all values are identical
		if (sd(df$Concentration) == 0) {
			return(data.frame(
				Test = "Constant data",
				P.value = 1,
				Significant.Pairs = "All values identical",
				stringsAsFactors = FALSE
			))
		}
		
		# Check normality only if n >= 3 and n <= 5000
		sample_size <- nrow(df)
		if (sample_size >= 3 && sample_size <= 5000) {
			norm_test <- shapiro.test(df$Concentration)
			norm_p <- norm_test$p.value
			normal <- norm_p > 0.05
		} else {
			normal <- FALSE
		}
		
		# Check homogeneity of variances
		var_test <- tryCatch({
			leveneTest(Concentration ~ Season, data = df)
		}, error = function(e) {
			return(list(`Pr(>F)` = c(NA, NA)))
		})
		
		var_p <- var_test$`Pr(>F)`[1]
		homogeneous <- ifelse(is.na(var_p), FALSE, var_p > 0.05)
		
		# Perform appropriate test
		if (normal && homogeneous) {
			# ANOVA
			anova_res <- aov(Concentration ~ Season, data = df)
			p_value <- summary(anova_res)[[1]]$`Pr(>F)`[1]
			test_used <- "ANOVA"
		} else {
			# Kruskal-Wallis
			kruskal_res <- kruskal.test(Concentration ~ Season, data = df)
			p_value <- kruskal_res$p.value
			test_used <- "Kruskal-Wallis"
		}
		
		# Initialize pairwise results
		pairwise_results <- "None"
		
		# Perform pairwise comparisons if overall significant
		if (!is.na(p_value) && p_value < 0.05) {
			if (test_used == "ANOVA") {
				pairwise <- tryCatch({
					TukeyHSD(anova_res)
				}, error = function(e) NULL)
				
				if (!is.null(pairwise)) {
					pairwise_df <- as.data.frame(pairwise$Season)
					sig_pairs <- pairwise_df %>%
						filter(`p adj` < 0.05) %>%
						rownames()
					pairwise_results <- if (length(sig_pairs) > 0) 
						paste(sig_pairs, collapse = "; ") else "No significant pairs"
				} else {
					pairwise_results <- "Tukey test failed"
				}
			} else {
				pairwise <- tryCatch({
					df %>% pairwise_wilcox_test(Concentration ~ Season, p.adjust.method = "bonferroni")
				}, error = function(e) NULL)
				
				if (!is.null(pairwise)) {
					sig_pairs <- pairwise %>%
						filter(p.adj < 0.05) %>%
						mutate(pair = paste0(group1, "-", group2)) %>%
						pull(pair)
					pairwise_results <- if (length(sig_pairs) > 0) 
						paste(sig_pairs, collapse = "; ") else "No significant pairs"
				} else {
					pairwise_results <- "Wilcox test failed"
				}
			}
		}
		
		# Return results
		data.frame(
			Test = test_used,
			P.value = p_value,
			Significant.Pairs = pairwise_results,
			stringsAsFactors = FALSE
		)
	}
	
	# Analyze data by Site and Metal
	results <- metal_long %>%
		group_by(Site, Metal) %>%
		group_modify(~ test_seasonal_differences(.x)) %>%
		ungroup()
	
	# Format and filter results
	final_results <- results %>%
		mutate(
			P.value = ifelse(Test %in% c("ANOVA", "Kruskal-Wallis"), round(P.value, 4), P.value),
			Result = case_when(
				Test == "Constant data" ~ "No variation (all values identical)",
				Test == "Insufficient data" ~ Significant.Pairs,
				P.value < 0.05 ~ "Significant seasonal difference",
				TRUE ~ "No significant difference"
			)
		) %>%
		select(Site, Metal, Test, P.value, Significant.Pairs, Result)
	
	# View significant results
	significant_results <- final_results %>%
		filter(Result == "Significant seasonal difference")
	
	# Print results
	print(significant_results, n = Inf)
	
	# export
	write_csv(
		final_results,
		"outputs/tables/final_results.csv"
	)
}

write_csv(
	metal_long,
	"outputs/tables/metal_long.csv"
)

{ # Are there significant metal differences between sites?
	
	# ---- helper: safe pairwise (base stats only) ----
	pairwise_safe <- function(df, method = "bonferroni") {
		ss <- df %>%
			dplyr::group_by(Site) %>%
			dplyr::summarise(n = dplyr::n(),
								  v = var(Concentration),
								  .groups = "drop")
		if (any(ss$n < 2L | ss$v == 0 | is.na(ss$v)))
			return("Skipped (constant or n<2)")
		
		out <- stats::pairwise.wilcox.test(df$Concentration,
													  df$Site,
													  p.adjust.method = method)
		sig <- which(out$p.value < 0.05, arr.ind = TRUE)
		if (!nrow(sig)) return("No significant pairs")
		
		apply(sig, 1, \(ix) paste(rownames(out$p.value)[ix[1]],
										  colnames(out$p.value)[ix[2]],
										  sep = "-")) |>
			paste(collapse = "; ")
	}
	
	# ---- helper: main test ----
	test_site_differences <- function(df) {
		
		df <- dplyr::filter(df, !is.na(Concentration))
		
		# need ≥2 sites
		if (dplyr::n_distinct(df$Site) < 2L)
			return(data.frame(Test = "Insufficient data",
									P.value = NA,
									Significant.Pairs = "≤1 site"))
		
		# flat slice
		if (length(unique(df$Concentration)) == 1L)
			return(data.frame(Test = "Constant data",
									P.value = 1,
									Significant.Pairs = "All values identical"))
		
		# check per-site n and variance
		ss <- df %>%
			dplyr::group_by(Site) %>%
			dplyr::summarise(n = dplyr::n(),
								  v = var(Concentration),
								  .groups = "drop")
		
		can_anova <- all(ss$n >= 2L) && all(ss$v > 0, na.rm = TRUE) &&
			shapiro.test(df$Concentration)$p.value > 0.05
		
		if (can_anova) {
			mod  <- oneway.test(Concentration ~ Site, data = df, var.equal = FALSE) # Welch
			pval <- mod$p.value
			test <- "Welch-ANOVA"
		} else {
			mod  <- kruskal.test(Concentration ~ Site, data = df)
			pval <- mod$p.value
			test <- "Kruskal-Wallis"
		}
		
		# pairwise only if overall sig
		pw <- "None"
		if (!is.na(pval) && pval < 0.05) {
			if (test == "Welch-ANOVA") {
				tuk <- TukeyHSD(aov(Concentration ~ Site, data = df))$Site
				sig <- rownames(tuk)[tuk[, "p adj"] < 0.05]
				pw  <- if (length(sig)) paste(sig, collapse = "; ") else "No significant pairs"
			} else {
				pw <- pairwise_safe(df)
			}
		}
		
		data.frame(Test = test,
					  P.value = round(pval, 4),
					  Significant.Pairs = pw)
	}
	
	# ---- run ----
	site_results <- metal_long %>%
		dplyr::group_by(Season, Metal) %>%
		dplyr::group_modify(~ test_site_differences(.x)) %>%
		dplyr::ungroup()
	
	# (optional) write_csv(site_results, "outputs/tables/site_results.csv")
	write_csv(
		site_results,
		"outputs/tables/site_results.csv"
	)
}

{ # GAM copper concentration forecast demo for Zoe ----
	copper16 <- metal_wide %>%
		filter(Site == 16) %>%
		select(
			Catchment,
			Site,
			year,
			Season,
			Copper
		)
	write_csv(
		copper16,
		"outputs/tables/copper16.csv"
	)
	
	season_levels <- c("Winter", "Spring", "Summer", "Autumn")
	copper16 <- copper16 %>% 
		arrange(year, factor(Season, levels = season_levels)) %>% 
		mutate(season_idx = row_number())      # 1, 2, 3 …
	
	# 2. fit GAM
	m <- gam(Copper ~ s(season_idx, k = 20, bs = "cs"),
				data = copper16,
				method = "REML")
	
	# 3. predict eight more seasons
	h <- 8
	future_idx <- max(copper16$season_idx) + seq_len(h)
	
	pred <- predict(m,
						 newdata = data.frame(season_idx = future_idx),
						 se.fit  = TRUE)
	
	future <- tibble(
		season_idx = future_idx,
		Copper     = NA_real_,
		fit        = pred$fit,
		lower      = pred$fit - 1.96 * pred$se.fit,
		upper      = pred$fit + 1.96 * pred$se.fit
	)
	
	# 4. plot
	ggplot() +
		geom_line(data = copper16,
					 aes(season_idx, Copper),
					 colour = "black",
					 linewidth = 0.4) +
		geom_ribbon(data = future,
						aes(season_idx, ymin = lower, ymax = upper),
						fill = "grey70",
						alpha = 0.4) +
		geom_line(data = future,
					 aes(season_idx, fit),
					 colour = "blue",
					 linewidth = 0.7,
					 linetype = "dashed") +
		labs(x = "Season index",
			  y = "Copper concentration",
			  title = "Observed and predicted copper at site 16") +
		theme_classic()
	
	## copper16 now contains historical rows with NA forecasts
	## and h extra rows with forecast values
}