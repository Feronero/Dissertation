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
