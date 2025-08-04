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
	library(BiocManager)
	library(ALDEx2)
}

{ # Sub-Sampled cross-catch PERMANOVAs (with nested Site term) ----
	# set taxonomic level
		i <- 4
	# min sites per catchment
		min_sites <- bact_data[[i]] %>%
			distinct(Catchment, Site) %>%
			count(Catchment) %>%
			pull(n) %>% min()
	# min samps per site
		min_per_site <- bact_data[[i]] %>%
			count(Catchment, Site) %>%
			group_by(Catchment) %>%
			summarise(m = min(n)) %>%
			pull(m) %>% min()
	# function for one balanced iteration
		one_run <- function() {
	# randomly subsample sites from each catchment
			samp_sites <- bact_data[[i]] %>%
				distinct(Catchment, Site) %>%
				group_by(Catchment) %>%
				slice_sample(n = min_sites) %>%
				ungroup()
	# randomly subsample min samples from each selected site
			samp <- bact_data[[i]] %>%
				semi_join(samp_sites, by = c("Catchment", "Site")) %>%
				group_by(Catchment, Site) %>%
				slice_sample(n = min_per_site) %>%
				ungroup()
	# establish permutation restrictions
			perm <- how(
				blocks = samp$Catchment,
				plots = Plots(strata = samp$Site, type = "free"),
				within = Within(type = "free"),
				nperm = 999
			)
	# run one PERMANOVA iteration
			mod <- adonis2(
				samp %>% select_after("invsimpson") ~
					Catchment + Catchment/Site + Site/Season +
					Temp + pH,
				data = samp,
				permutations = perm,
				method = "bray",
				by = "terms",
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
		summary_stats <- out %>%
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
summary_stats
saveRDS(summary_stats, "outputs/backups/summary_stats.rds")
summary_stats <- readRDS("outputs/backups/summary_stats.rds")

{ # Sub-Sampled cross-catch PERMANOVAs (with rivergradient and veg surround) ----
	# set taxonomic level
		i <- 4
	# min sites per catchment
		min_sites <- bact_data[[i]] %>%
			distinct(Catchment, Site) %>%
			count(Catchment) %>%
			pull(n) %>% min()
	# min samps per site
		min_per_site <- bact_data[[i]] %>%
			count(Catchment, Site) %>%
			group_by(Catchment) %>%
			summarise(m = min(n)) %>%
			pull(m) %>% min()
	# function for one balanced iteration
		one_run <- function() {
	# randomly subsample sites from each catchment
			samp_sites <- bact_data[[i]] %>%
				distinct(Catchment, Site) %>%
				group_by(Catchment) %>%
				slice_sample(n = min_sites) %>%
				ungroup()
	# randomly subsample min samples from each selected site
			samp <- bact_data[[i]] %>%
				semi_join(samp_sites, by = c("Catchment", "Site")) %>%
				group_by(Catchment, Site) %>%
				slice_sample(n = min_per_site) %>%
				ungroup()
	# establish permutation restrictions
			perm <- how(
				blocks = samp$Catchment,
				plots = Plots(strata = samp$Site, type = "free"),
				within = Within(type = "free"),
				nperm = 999
			)
	# run one PERMANOVA iteration
			mod <- adonis2(
				samp %>% select_after("invsimpson") ~
					Catchment + Site + River_Gradient + Veg_Surround + Season +
					Temp + pH,
				data = samp,s
				permutations = perm,
				method = "bray",
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
		catch_non_nested <- out %>%
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
catch_non_nested

{ # Sub-Sampled Catchment-level PERMANOVAs ----
	# results storage lists
		percatch_raw <- list()
		percatch_summary <- list()
	# loop over each catchment
		for (j in keyvals$catchments) {
	# set taxonomic level
			i <- 4
	# create subset for current catchment
			catch_subset <- bact_data[[i]] %>%
				filter(Catchment == j)
	# min samps per site
			min_per_site <- catch_subset %>%
				count(Site) %>%
				pull(n) %>%
				min()
	# function for one balanced iteration
			one_run <- function(iter) {
	# randomly subsample min samples from each selected site
				samp <- catch_subset %>%
					group_by(Site) %>%
					slice_sample(n = min_per_site) %>%
					ungroup()
	# establish permutation restrictions
				perm <- how(
					plots = Plots(strata = samp$Site, type = "free"),
					within = Within(type = "free"),
					nperm = 999
				)
	# run one PERMANOVA iteration
				mod <- adonis2(
					samp %>% select_after("invsimpson") ~
						Site + Site/Season +
						Temp + pH,
					data = samp,
					permutations = perm,
					method = "bray",
					by = "terms",
					parallel = 2
				)
	# accumulate key output from each iteration
				tibble(
					Catchment = j,
					Iteration = iter,
					Term = rownames(mod),
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
				~ one_run(.x),
				.options = furrr_options(seed = TRUE),
				.progress = TRUE
			)
	# save raw results of all runs
			percatch_raw[[j]] <- out
		}
		percatch_summary <- percatch_raw %>%
			bind_rows() %>%
			filter(!is.na(p)) %>%
			group_by(Catchment, Term) %>%
			summarise(
				Iterations = max(Iteration),
				lo_R2 = quantile(R2, .025),
				hi_R2 = quantile(R2, .975),
				median_R2 = median(R2),
				median_p = median(p),
				`%p<0.1` = mean(p < 0.1),
				`%p<0.05` = mean(p < 0.05) * 100,
				`%p<0.01` = mean(p < 0.01) * 100,
				`%p<0.001` = mean(p < 0.001) * 100,
				.groups = "drop"
			) %>%
			group_by(Catchment) %>%
			mutate(p_holm = p.adjust(median_p, method = "holm")) %>%
			split(.$Catchment)
}
percatch_summary
saveRDS(percatch_summary, "outputs/backups/percatch_summary.rds")
percatch_summary <- readRDS("outputs/backups/percatch_summary.rds")

{ # Pairwise PERMANOVA within each catchment ----
		set.seed(42)
		plan(multisession, workers = 8)
	# builds a tibble listing all site-pairs within all catchments
		pairs_df <- keyvals$catchments %>% 
		  map_dfr(function(i) {
		    sites <- bact_data[[4]] %>%
		    	filter(Catchment == i) %>%
		    	pull(Site) %>%
		    	unique() %>%
		    	as.numeric()
		    combn(sites, 2, simplify = FALSE) %>%
		    	map_df(~{
		      	these <- .x
		      	min_samps <- bact_data[[4]] %>% 
		      		filter(Site %in% these) %>%
		      		count(Site) %>%
		      		pull(n) %>%
		      		min()
		      	tibble(
		      		Catchment = i,
		            Pair = paste(these, collapse = ","),
		            Site1 = these[1],
		            Site2 = these[2],
		            min_samps = min_samps
		      	)
		    })
		  })
	# function that performs all iterations for one pair
		one_pair_runs <- function(row) {
			i <- row$Catchment
			s1 <- row$Site1
			s2 <- row$Site2
			pair_label <- row$Pair
			min_samps <- row$min_samps
			map_dfr(1:50, function(iter) {
				pair_data <- bact_data[[4]] %>% 
					filter(Site %in% c(s1, s2)) %>% 
					group_by(Site) %>%
					slice_sample(n = min_samps) %>%
					ungroup()
				mod <- adonis2(
					pair_data %>% select_after("invsimpson") ~
						Site,
					method = "bray",
					by = "terms",
					data = pair_data,
					permutations = 999
				)
				tibble(
					Catchment = i,
					Pair = pair_label,
					Site1 = s1,
					Site2 = s2,
					Iteration = iter,
					R2 = mod$R2[1],
					p = mod$`Pr(>F)`[1]
				)
			})
		}
	# run in parallel, each worker is assigned to one pair and handles each iteration
		raw_results <- future_map_dfr(
			split(pairs_df, 1:nrow(pairs_df)),
			~ one_pair_runs(.x),
			.options = furrr_options(seed = TRUE),
			.progress = TRUE
		)
	# summarise the raw results
		persite_summary <- raw_results %>%
			group_by(Catchment, Pair) %>%
			summarise(
				Iterations = max(Iteration),
				Site1 = first(Site1),
				Site2 = first(Site2),
				lo_R2 = quantile(R2, .025),
				hi_R2 = quantile(R2, .975),
				median_R2 = median(R2),
				median_p = median(p),
				`%p<0.1` = mean(p < 0.1) * 100,
				`%p<0.05` = mean(p < 0.05) * 100,
				`%p<0.01` = mean(p < 0.01) * 100,
				`%p<0.001` = mean(p < 0.001) * 100,
				.groups = "drop"
			) %>%
			group_by(Catchment) %>%
			mutate(p_holm = p.adjust(median_p, method = "holm")) %>%
			ungroup()
	# drop non-sig results
		persite_valid <- persite_summary %>%
			filter(p_holm < 0.1, `%p<0.05` > 80)
		persite_valid_splitbycatch <- persite_valid %>%
			group_by(Catchment) %>%
			split(.$Catchment)
}
persite_valid
saveRDS(persite_valid, "outputs/backups/persite_valid.rds")
persite_valid <- readRDS("outputs/backups/persite_valid.rds")

persite_valid_splitbycatch
saveRDS(persite_valid_splitbycatch, "outputs/backups/persite_valid_splitbycatch.rds")
persite_valid_splitbycatch <- readRDS("outputs/backups/persite_valid_splitbycatch.rds")

{ # ALDEx2 in parallel via future ----
results_list <- list()
bact_matrix <- bact_data[[4]] %>%
	column_to_rownames(var = "Sample") %>%
	select_after("invsimpson") %>%
	t()
bact_meta <- bact_data[[4]] %>%
	column_to_rownames(var = "Sample") %>%
	select(Catchment, Site) %>%
	as.data.frame()
pairwise_signif <- persite_valid %>%
	select(Site1, Site2)
	set.seed(48)
	plan(multisession, workers = 8)
		process_pair <- function(s1, s2, bact_matrix, bact_meta) {
			sel <- bact_meta$Site %in% c(s1, s2)
			samples <- rownames(bact_meta)[sel]
			otu <- bact_matrix[, samples, drop = FALSE]
			groups <- as.numeric(droplevels(bact_meta$Site[sel]))
			otu <- otu[rowSums(otu) > 0, , drop = FALSE]
			x.clr    <- aldex.clr(otu, conds = groups, mc.samples = 128, denom = "iqlr", verbose = FALSE)
			x.kw     <- aldex.kw(x.clr, verbose = TRUE)
			x.effect <- aldex.effect(x.clr, CI = TRUE, verbose = TRUE)
			out <- data.frame(
				feature  = rownames(x.kw),
				kws.p    = x.kw$kw.ep,
				kws.q    = x.kw$kw.eBH,
				effect   = x.effect$effect,
				diff.btw = x.effect$diff.btw,
				overlap  = x.effect$overlap,
				stringsAsFactors = FALSE
			)
			out[order(-abs(out$effect), out$kws.q), ]
		}
ALDEx2_results <- future_map2(
  pairwise_signif$Site1,
  pairwise_signif$Site2,
  process_pair,
  bact_matrix = bact_matrix,
  bact_meta   = bact_meta,
  .options    = furrr_options(seed = TRUE),
  .progress = TRUE
) %>%
  set_names(paste0(pairwise_signif$Site1, "_vs_", pairwise_signif$Site2))
}
View(ALDEx2_results)
saveRDS(ALDEx2_results, "outputs/backups/ALDEx2_results.rds")
ALDEx2_results <- readRDS("outputs/backups/ALDEx2_results.rds")


df_all <- enframe(results_list, name = "comparison", value = "tbl") %>%
  unnest_longer(tbl) %>%       # this gives one row per feature per comparison
  unnest_wider(tbl)

df_sig <- df_all %>%
  filter(
  	kws.q < 0.05,
	abs(effect) > 1,
   overlap < 0.25) %>%
	arrange(comparison, desc(abs(effect)))
