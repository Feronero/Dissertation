{ # Compare PERMANOVA performance with + without Temp ----
	# Fit both PERMANOVA models
perm <- how(
	blocks = bact_data[[4]]$Catchment,
	plots = Plots(
		strata = bact_data[[4]]$Site,
		type = "free"
	),
	within = Within(type = "free"),
	nperm = 999
)

mod_with <- adonis2(
	bact_data[[4]] %>% select_after("invsimpson") ~
		Catchment + Catchment/Site + Site/Season + Veg_Surround + pH + Temp,
	data = bact_data[[4]],
	permutations = perm,
	by = "margin"
)
	
		mod_with <- adonis2(
			bact_data[[4]] %>% select_after("invsimpson") ~
				Catchment / Site + Catchent / Season + River_Gradient + Veg_Surround + pH + Temp,
			data = bact_data[[4]],
			method = "bray",
			by = "margin",
			strata = structure
		)
		mod_no <- adonis2(
			bact_data[[4]] %>% select_after("invsimpson") ~
				Catchment + Season + River_Gradient + Veg_Surround + pH,
			data = bact_data[[4]],
			method = "bray",
			by = "margin",
			strata = bact_data[[4]]$Batch
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

{ # Sub-sampled cross-catch PERMANOVA ----
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
			1:20,
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

{ # Test for sig. diff in dispersion between groups ----
one_disp_run <- function() {
  # balanced subsampling as before
  samp_sites <- bact_data[[i]] %>%
    distinct(Catchment, Site) %>%
    group_by(Catchment) %>%
    slice_sample(n = min_sites) %>%
    ungroup()
  samp <- bact_data[[i]] %>%
    semi_join(samp_sites, by = c("Catchment","Site")) %>%
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

  dist_mat <- vegdist(samp %>% select_after("invsimpson"), method = "bray")

  out <- bind_rows(

    # 1. Catchment-level dispersion
	{
      bd_c <- betadisper(dist_mat, samp$Catchment)
      pt_c <- permutest(bd_c, permutations = perm)
      tibble(factor="Catchment", level="Catchment", F.obs=pt_c$tab$F[1], p.value=pt_c$tab$`Pr(>F)`[1])
   },

    # 2. Site-level dispersion (sites nested within catchments)
    {
      revo <- samp %>% mutate(Catch_Site = interaction(Catchment, Site, sep="_"))
      bd_s <- betadisper(dist_mat, revo$Catch_Site)
      pt_s <- permutest(bd_s, permutations = perm)
      tibble(factor="Site", level="Site_within_Catchment", F.obs=pt_s$tab$F[1], p.value=pt_s$tab$`Pr(>F)`[1])
    },

    # 3. Season-level dispersion (season nested within Site)
    {
      revo <- samp %>%
        mutate(Catch_Site_Season = interaction(Catchment, Site, Season, sep="_"))
      bd_se <- betadisper(dist_mat, revo$Catch_Site_Season)
      pt_se <- permutest(bd_se, permutations = perm)
      tibble(factor="Season", level="Season_within_Site", F.obs=pt_se$tab$F[1], p.value=pt_se$tab$`Pr(>F)`[1])
    }
  )

  return(out)
}
	
	set.seed(42)
plan(multisession, workers=8)
disp_out <- future_map_dfr(
	1:8,
	~ one_disp_run(),
	.options = furrr_options(seed = TRUE),
	.progress = TRUE
)

disp_summary <- disp_out %>%
 group_by(factor, level) %>%
 summarise(
   median_F = median(F.obs),
   lo_F = quantile(F.obs, 0.025),
   hi_F = quantile(F.obs, 0.975),
   `%p<0.05` = mean(p.value < 0.05)*100,
   `%p<0.01` = mean(p.value < 0.01)*100,
   `%p<0.001` = mean(p.value < 0.001)*100,
   .groups="drop"
 )
}

{ # View Final PERMANOVA results ----
	summary_stats
}

{ # Catchment-level PERMANOVAs ----
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
		
		percatch_adj <- percatch_raw %>%
			bind_rows() %>%
			filter(!is.na(p)) %>%
			group_by(Catchment) %>%
			mutate(p_holm = p.adjust(p, method = "holm")) %>%
			ungroup()
		
		for (j in keyvals$catchments) {
			percatch_summary[[j]] <- percatch_adj %>%
				filter(Catchment == j) %>%
				group_by(Term) %>%
				summarise(
					median_R2 = median(R2),
					lo_R2 = quantile(R2, .025),
					hi_R2 = quantile(R2, .975),
					`%p_holm<0.05` = median(p_holm < 0.05) * 100,
					`%p_holm<0.01` = median(p_holm < 0.01) * 100,
					`%p_holm<0.001` = median(p_holm < 0.001) * 100,
					.groups = "drop"
				)
		}
}

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
		summarised_results <- raw_results %>%
			group_by(Catchment, Pair) %>%
			summarise(
				Site1 = first(Site1),
				Site2 = first(Site2),
				median_p = median(p),
				prop_sig = mean(p < 0.05),
				.groups = "drop"
			)
	# apply Holm correction by catchment
		corrected_results <- summarised_results %>%
			group_by(Catchment) %>%
			mutate(p_holm = p.adjust(median_p, method = "holm")) %>%
			ungroup()
	# drop non-sig results
		valid_results <- corrected_results %>%
			filter(p_holm < 0.1 & prop_sig >= 0.8)
}

{ # ANCOMB2 Attempt. Terrible. ----
bact_matrix <- bact_data[[4]] %>%
	column_to_rownames(var = "Sample") %>%
	select_after("invsimpson") %>%
		as.data.frame()
bact_meta <- bact_data[[4]] %>%
	column_to_rownames(var = "Sample") %>%
	select(Catchment, Site) %>%
	as.data.frame()
	
output <- ancombc2(
	data = bact_matrix,
	taxa_are_rows = FALSE,
	meta_data = bact_meta,
	fix_formula = "Site",
	rand_formula = "(1 | Catchment / Site)",
	p_adj_method = "holm",
	group = "Site",
	prv_cut = 0.10,
	lib_cut = 0,
	pseudo_sens = TRUE,
	struc_zero = TRUE,
	neg_lb = TRUE,
	global = TRUE,
	pairwise = TRUE,
	alpha = 0.05,
	n_cl = 8
)
res_primary   <- output$res_prim
res_global    <- output$res_global
res_pairwise  <- output$res_pair

 # attempt to fix non-convergence
	# example: select the top 10 most abundant taxa
taxa_to_test <- names(sort(colSums(bact_matrix), decreasing = TRUE))[1:10]

pilot_counts <- bact_matrix[, taxa_to_test]
pilot_meta   <- bact_meta

pseudo = 1
pilot_df <- pilot_counts %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(cols = -Sample, names_to = "Taxon", values_to = "count") %>%
  mutate(Y = log(count + pseudo)) %>%
  left_join(bact_meta %>% rownames_to_column("Sample"), by = "Sample") %>%
  mutate(Site = factor(Site), Catchment = factor(Catchment))

library(lme4)
pilot_mod <- lmer(Y ~ Site + (1 | Catchment / Site), data = pilot_df, REML = TRUE)

check_model(pilot_mod)
test_performance(pilot_mod)
simulateResiduals(pilot_mod, plot = TRUE)

is_singular <- isSingular(pilot_mod, tol = 1e-4)
print(is_singular)

library(performance)
check_singularity(pilot_mod)  # term-level singularity report

summary(pilot_mod)

all_fits <- allFit(pilot_mod)
summary(all_fits)

ss <- summary(all_fits)
ss$fixef
ss$sdcor
ss$llik
}

{ # ALDEx2 in series ----
results_list <- list()
	
bact_matrix <- bact_data[[4]] %>%
	column_to_rownames(var = "Sample") %>%
	select_after("invsimpson") %>%
	t()
bact_meta <- bact_data[[4]] %>%
	column_to_rownames(var = "Sample") %>%
	select(Catchment, Site) %>%
	as.data.frame()
pairwise_signif <- valid_results %>%
	select(Site1, Site2)

for (i in 1:nrow(pairwise_signif)) {
  s1 <- pairwise_signif$Site1[i]
  s2 <- pairwise_signif$Site2[i]
  
  # Identify samples in these two site groups
  sel <- bact_meta$Site %in% c(s1, s2)
  samples <- rownames(bact_meta)[sel]
  otu <- bact_matrix[, samples, drop = FALSE]
  groups <- bact_meta$Site[sel] %>%
  	droplevels() %>%
  	as.numeric()
  
  # You may want to filter features with zero-sum or low prevalence
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]
  
  set.seed(123)
  # Step 1: CLR with Monte Carlo sampling
  x.clr <- aldex.clr(otu, conds = groups,
                     mc.samples = 128,
                     denom = "iqlr",
                     verbose = FALSE)
  
  # Step 2: Kruskal-Wallis test across groups
  x.kw <- aldex.kw(x.clr, verbose = TRUE)
  
  # Step 3: calculate effect sizes (with confidence intervals if desired)
  x.effect <- aldex.effect(x.clr, CI = TRUE, verbose = TRUE)
  
  # Combine results
out <- data.frame(
	feature = rownames(x.kw),
	kws.p = x.kw$kw.ep,
	kws.q = x.kw$kw.eBH,
	effect = x.effect$effect,
	diff.btw = x.effect$diff.btw,
	overlap = x.effect$overlap,
	stringsAsFactors = FALSE
)
  # Order by effect size and FDR
	  out <- out[ order(-abs(out$effect), out$kws.q), ]
  
  results_list[[paste0(s1, "_vs_", s2)]] <- out
}

# Access results for a given pair:
results_list[["SiteA_vs_SiteB"]]
}

{ # ALDEx2 in parallel ----
	plan(multisession, workers = availableCores())
		process_pair <- function(s1, s2, bact_matrix, bact_meta) {
			sel <- bact_meta$Site %in% c(s1, s2)
			samples <- rownames(bact_meta)[sel]
			otu <- bact_matrix[, samples, drop = FALSE]
			groups <- as.numeric(droplevels(bact_meta$Site[sel]))
			otu <- otu[rowSums(otu) > 0, , drop = FALSE]
			set.seed(123)
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
results_list <- future_map2(
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
		comm <- bact_data[[4]] %>% select_after("invsimpson")
	# copy data
		data <- bact_data[[4]]
	# Store db-RDA plots
		db_rda_plots <- list()
	# Factors to test
		factors <- c("Catchment", "Season", "River_Gradient", "Veg_Surround")
	# Loop through each factor
		for (var in factors) {
	# If numeric, bin to factor for plotting
		  var_data <- data[[var]]
		  if (is.numeric(var_data)) {
		    data[[paste0(var, "_binned")]] <- cut(var_data, breaks = 4)
		    var_plot <- paste0(var, "_binned")
		  } else {
		    var_plot <- var
		  }
	# Model
		  form <- as.formula(paste0("comm ~ ", var))
		  mod <- capscale(form, data = data, distance = "bray")
	# Scores
		  scores_df <- data %>%
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
		      aes(colour = group, fill = group, shape = group),
		      size = 3,
		      alpha = 0.8
		   ) +
			scale_shape_manual(values = 21:24) +
		    theme_classic() +
		    labs(
		      x = "db-RDA Axis 1",
		      y = "db-RDA Axis 2",
		      title = paste("db-RDA:", var)
		    ) +
		    guides()  # hide fill legend if it duplicates colour
		  
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
		        geom = "tile", contour = FALSE, n = 800,
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

{ # per-catchment PERMANOVAs ----
	
	percatch_perm <- list()
	
	for (i in keyvals$catchments) {
		
		catch_subset <- bact_data[[4]] %>%
			filter(Catchment == i)
		catch_bact <- catch_subset %>%
			select_after("invsimpson")
		
		percatch_perm[[i]] <- adonis2(
			catch_bact ~
				Season + River_Gradient + pH + Veg_Surround,
			data = catch_subset,
			strata = catch_subset$Batch,
			by = "margin"
		)
		percatch_perm[[i]] %>%
			print()
	}
}

{ # Attempt to fit parsimonious db-RDAs unique to each river ----
		selection <- list()
		for (i in keyvals$catchments) {
			catch_subset <- filter(bact_data[[4]], Catchment == i)
			comm <- catch_subset %>% select_after("invsimpson")
			env <- catch_subset %>% select(Season, River_Gradient, Veg_Surround, pH)
			
			null_mod <- capscale(comm ~ 1, data = env, distance = "bray")
			full_mod <- capscale(comm ~ Season + River_Gradient + Veg_Surround + pH,
			                    data = env, distance = "bray")
			
			opt <- ordiR2step(null_mod, scope = formula(full_mod), R2scope = TRUE,
			  					direction = "both",
			                 permutations = how(nperm = 499), trace = FALSE)
			
			selection[[i]] <- list(model = opt,
			                      anova_terms = anova(opt, by = "terms", permutations = 999))
			results[[i]] <- data.frame(
				Catchment = i,
				Formula = formula(selection[[i]]$model) %>%
					deparse(),
				n_samps = nrow(catch_subset)
			)
		}
}

formula(results[[1]]$model) %>%
	deparse()

{ # misc ----
for (i in keyvals$catchments) {
	bact_collapsed[[4]] %>%
		filter(Catchment == i) %>%
		nrow() %>%
		print()
}
}
