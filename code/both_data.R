{ # Merge metal PCs with bacteria data ----
	master <- list()
	for (i in 1:keyvals$n_levels) {
		master[[i]] <- merge(
			x = bact_data[[i]],
			y = meta_pc_agg,
			by = c("Site", "Season")
		) %>%
			relocate(
				starts_with("PC"),
				.after = "invsimpson"
			)
	}
	master_collapsed <- list()
	for (i in 1:keyvals$n_levels) {
		master_collapsed[[i]] <- merge(
			x = bact_collapsed[[i]],
			y = meta_pc_agg,
			by = c("Site", "Season")
		) %>%
			relocate(
				starts_with("PC"),
				.after = "invsimpson"
			)
	}
	write_csv(
		master[[1]],
		"outputs/tables/master[[1]].csv"
	)
	write_csv(
		master_collapsed[[1]],
		"outputs/tables/master_collapsed[[1]].csv"
	)
}

{ # Drop Hayle and Cober - only 7 samples in master_collapsed each, unable to safely run rarefied subsampling verification runs ----
	CarnonRed_collapsed <- list()
	for (i in 1:keyvals$n_levels) {
		CarnonRed_collapsed[[i]] <- master_collapsed[[i]] %>% filter(!Catchment %in% c("Hayle", "Cober"))
	}
	CarnonRed_OTUs <- CarnonRed_collapsed[[4]] %>% select_after("PC3")
}

{ # Test PERMANOVA variables for overdispersion ----
	# full bact data, no metals
	
		# define function to test dispersion change across a cont. variable
		perm_test <- function(y, x, nperm=999) {
			f0 <- summary(lm(y ~ x))$fstatistic[1]
			permF <- replicate(nperm, summary(lm(sample(y) ~ x))$fstatistic[1])
			mean(permF >= f0)
		}	
	
		# separate bacterial data from metadata, calc distances
			bact_OTUs <- select_after(bact_collapsed[[1]], "invsimpson")
			bact_dist <- bact_OTUs %>%
				vegdist(method = "bray")
		# test each categorical factor for het. dispersion
			permutest(betadisper(bact_dist, group = bact_collapsed[[4]]$Catchment), permutations = 999)
			permutest(betadisper(bact_dist, group = bact_collapsed[[4]]$Site), permutations = 999)
			permutest(betadisper(bact_dist, group = bact_collapsed[[4]]$Season), permutations = 999)
			permutest(betadisper(bact_dist, group = bact_collapsed[[4]]$Veg_Surround), permutations = 999)
		# test each continuous factor for het. dispersion
			# run betadisper on all samples using a dummy group
				distances <- betadisper(bact_dist, group = rep(1, nrow(bact_OTUs)))$distances
			# calculate results for each continuous variable
				perm_test(distances, bact_collapsed[[4]]$River_Gradient)
				perm_test(distances, bact_collapsed[[4]]$Temp)
				perm_test(distances, bact_collapsed[[4]]$pH)
	# bact + metal data
		# separate bacterial data from metadata, calc distances
			master_OTUs <- select_after(master_collapsed[[1]], "PC3")
			master_dist <- master_OTUs %>%
				vegdist(method = "bray")
		# test each categorical factor for het. dispersion
			permutest(betadisper(master_dist, group = master_collapsed[[4]]$Catchment), permutations = 999)
			permutest(betadisper(master_dist, group = master_collapsed[[4]]$Site), permutations = 999)
			permutest(betadisper(master_dist, group = master_collapsed[[4]]$Season), permutations = 999)
			permutest(betadisper(master_dist, group = master_collapsed[[4]]$Veg_Surround), permutations = 999)
		# test each continuous factor for het. dispersion
		# run betadisper on all samples using a dummy group
			distances <- betadisper(master_dist, group = rep(1, nrow(master_OTUs)))$distances
		# calculate results for each continuous variable
			perm_test(distances, master_collapsed[[4]]$River_Gradient)
			perm_test(distances, master_collapsed[[4]]$Temp)
			perm_test(distances, master_collapsed[[4]]$pH)
			
		# clean up
			rm(distances)
}

{ # Test PERMANOVA variables for VIF ---
	vif_vals <- vif(lm(
		Temp ~ pH + Site + Catchment + River_Gradient + Season,
		data = bact_collapsed[[4]]
	))
	print(vif_vals)
}

{ # Overall PERMANOVAs (bacteria only) ----
 
}

{ # Overall PERMANOVAs (metal + bact) ----
	# all metal+bact data (cannot be verified via subsampling; insufficient samples from Hayle + Cober)
		# Identify variables with high R2
			full_mod <- capscale(
			master_collapsed_OTUs ~
				Catchment + Season + River_Gradient + Veg_Surround + PC1 + PC2 + PC3,
			data = master_collapsed[[4]]
			)
			null_mod <- capscale(
				master_collapsed_OTUs ~ 1,
				data = master_collapsed[[4]]
			)
			
			sel <- ordiR2step(
				null_mod,
				scope = formula(full_mod),
				direction = "both",
				R2scope = TRUE,
				permutations = how(
					nperm = 999,
					blocks = master_collapsed[[4]]$Site
				)
			)
			sel
			summary(sel$anova)
		# run PERMANOVA
		# 57 samples; max 5 terms
			for (i in 4) {						# ready for multiple tax. levels if needed
				result <- adonis2(
					master_collapsed[[i]] %>% select_after("PC3") ~
						Catchment + Season + PC1 + PC2 + PC3,
					data = master_collapsed[[i]],
					method = "bray",
					by = "margin",
					strata = master_collapsed[[i]]$Site
				)
				print(result)
			}
	# CarnonRed subset (Hayle and Cober removed). Can be verified via subsampling
		# Identify variables with high R2
			full_mod <- capscale(
				CarnonRed_collapsed[[4]] %>% select_after("PC3") ~
					Season + River_Gradient + Veg_Surround + PC1 + PC2 + PC3,
				data = CarnonRed_collapsed[[4]]
			)
			null_mod <- capscale(
				CarnonRed_collapsed[[4]] %>% select_after("PC3") ~
					1,
				data = CarnonRed_collapsed[[4]]
			)
			
			sel <- ordiR2step(
				null_mod,
				scope = formula(full_mod),
				direction = "both",
				R2scope = TRUE,
				permutations = how(
					nperm = 999,
					blocks = CarnonRed_collapsed[[4]]$Site
				)
			)
			sel
			summary(sel$anova)		
		# Run PERMANOVA
		# 43 samples; max 4 terms
		# ordiR2step selected only Season
			for (i in 4) {						# ready for multiple tax. levels if needed
				result <- adonis2(
					CarnonRed_collapsed[[i]] %>% select_after("PC3") ~
						Season + PC1 + PC2 + PC3,
					data = CarnonRed_collapsed[[i]],
					method = "bray",
					by = "margin",
					strata = CarnonRed_collapsed[[i]]$Site
				)
				print(result)
			}
}

{ # overall alpha-div (bact only)
	cor.test(
		bact_collapsed[[4]]$invsimpson,
		bact_collapsed[[4]]$River_Gradient,
		method = "spearman"
	)
}

{ # cross-catchment db-RDA
# db‑RDA
mod_all <- capscale(
	vegdist(bact_collapsed[[4]] %>% select_after("invsimpson"), method = "bray") ~
		Catchment + Season + Veg_Surround + River_Gradient + Temp + pH,
	data = bact_collapsed[[4]]
)

# test
anova(
	mod_all,
	by = "terms",
	permutations = how(nperm = 999, blocks = bact_collapsed[[4]]$Site)
)

sig_cont_vars <- c("River_Gradient")
sig_cat_vars <- c("v")

# Site scores with metadata
site_scores <- scores(mod_all, display = "sites") %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  bind_cols(bact_collapsed[[4]] %>% select(Catchment, Season))

# Biplot (environmental variable) scores
biplot_scores <- scores(mod_all, display = "bp") %>%
	as.data.frame() %>%
	mutate(Variable = rownames(.)) %>%
	filter(Variable %in% sig_cont_vars)

# plot
ggplot() +
  geom_hline(yintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_point(data = site_scores, aes(CAP1, CAP2, color = Catchment, shape = Season), size = 3) +
  geom_segment(data = biplot_scores, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  geom_text(data = biplot_scores, aes(x = CAP1, y = CAP2, label = Variable),
            size = 3, hjust = 1.1, vjust = 1.1) +
  theme_minimal() +
  labs(title = "db-RDA of bacterial communities (All Rivers)",
       x = "db-RDA1", y = "db-RDA2") +
  theme(panel.grid = element_blank())
}

{ # db-RDA v2 ----
	# If you want a separate plot for each variable (site points remain same):
for (var in sig_cont_vars) {
  biplot_single <- biplot_scores %>% filter(Variable == var)
  
  ggplot() +
    geom_hline(yintercept = 0, color = "gray80") +
    geom_vline(xintercept = 0, color = "gray80") +
    geom_point(data = site_scores, aes(CAP1, CAP2, color = Catchment), size = 2) +
    geom_segment(data = biplot_single, aes(x = 0, y = 0, xend = CAP1, yend = CAP2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    geom_text(data = biplot_single, aes(x = CAP1, y = CAP2, label = Variable),
              hjust = 1.1, vjust = 1.1, size = 3) +
    theme_minimal() +
    labs(title = paste("db-RDA: Effect of", var),
         x = "db-RDA1", y = "db-RDA2") +
    theme(panel.grid = element_blank()) -> p
  print(p)
}
}






{ # Identify top phyla ----
	# select bacterial columns
	bact_cols <- names(bact_data[[1]][, 22:ncol(bact_data[[1]])])
	# converts raw sequence counts to relative abundance percentages
	df_agg <- bact_data[[1]] %>%
		mutate(across(
			all_of(bact_cols),
			~ .x / rowSums(across(all_of(bact_cols)), na.rm = TRUE) * 100
		)) %>%
		group_by(Catchment, Season) %>%
		summarise(across(all_of(bact_cols), \(x) mean(x, na.rm = TRUE))) %>%
		ungroup()
	
	top_phyla <- df_agg %>%
		summarise(across(-c(Catchment, Season), mean)) %>%
		pivot_longer(everything(), names_to = "phylum", values_to = "abundance") %>%
		arrange(desc(abundance)) %>%
		slice_head(n = 10) %>%
		pull(phylum)
	
	# Combine non-top phyla into "Other"
	df_agg_top <- df_agg %>%
		mutate(
			Other = rowSums(
				across(
					!all_of(
						c("Catchment", "Season", top_phyla)
					)
				)
			)
		) %>%
		select(
			all_of(
				c("Catchment", "Season", top_phyla, "Other")
			)
		) %>%
		pivot_longer(
			-c(Catchment, Season),
			names_to = "phylum",
			values_to = "abundance"
		)
	# Order seasons chronologically
#	season_order <- c("January", "April", "July", "October")
#	df_agg_top$Season <- factor(
#		x = df_agg_top$Season,
#		levels = season_order
#	)
}

{ # Bar plot of top phyla by catchment and season ----
	ggplot(
		df_agg_top, 
		aes(
			x = Season,  # Use Month as x-axis (no interaction with Catchment)
			y = abundance, 
			fill = phylum
		)
	) +
		geom_col(
			position = "stack",
			width = 0.75
		) +
		# separate rivers on separate facets
		facet_wrap(
			~ Catchment, 
			scales = "free_x",  # Let x-axis vary per river
			strip.position = "top",  # Place river labels at top
			nrow = 1  # Force all rivers into one row
		) +
		scale_fill_viridis_d(
			option = "D",
			labels = function(x) {
				# Create expression for each label
				sapply(
					x,
					function(phylum) {
						if (phylum == "Other") {
							expression(plain(Other))
						} else {
							# Extract just the phylum name (after last '__')
							phylum_name <- sub("^.*__", "", phylum)
							bquote(italic(.(phylum_name)))
						}
					},
					simplify = FALSE
				)
			}
		) +
		labs(
			x = "Season",  # Now x-axis shows months (rivers are in facets)
			y = "Relative Abundance (%)",
			fill = "Phylum",
			title = "Top 10 Bacterial Phyla by River and Season"
		) +
		theme_classic() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			#			axis.title.x = element_text(margin = margin(t = 10)),
			panel.spacing = unit(1, "lines"),  # Space between river clusters
			strip.background = element_blank(),  # Remove facet background
			strip.text = element_text(size = 10, face = "bold"),  # Style river labels
			legend.position = "bottom",
			plot.title = element_text(hjust = 0.5),
			legend.title = element_text(angle = 90, hjust = 0.5),
			legend.spacing.x = unit(50, "cm"),
			legend.box.margin = margin(0,0,0,0)
		)
}

{ # a-Diversity changes? ----
	# check richness distribution
	qqnorm(bact_data[[4]]$richness)
	qqline(bact_data[[4]]$richness)
	shapiro_test(bact_data[[4]]$richness)
	hist(bact_data[[4]]$richness)
}

{ # Abundance changes? ----
	
}

{ # Indicator species? ----
	ind_result <- multipatt(
		master[[1]][, 24:ncol(master[[1]])],
		master[[1]]$Catchment,
		control = how(nperm = 999)
	)
	ind_result <- ind_result$sign %>%
		filter(p.value <= 0.05)
	
	View(ind_result)
}

{ # Envfit ----
	comm <- df[, 24:ncol(df)]  # adjust col numbers if needed
	env <- df[, 1:23]          # adjust col numbers if needed
	
	# Run capscale
	mod <- capscale(
		master[[4]][, 24:ncol(master[[4]])] ~
			Catchment + Alt_Class + Veg_Surround + Season + Temp + pH + PC1 + PC2 + PC3,
		data = master[[4]],
		distance = "bray"
	)
	
	# Check the result
	print(mod)
	anova(mod, by = "term", permutations = 999)  # test significance of terms
	
	# Plot ordination
	plot(mod)
	
	# Fit env variables to ordination (optional)
	ef <- envfit(
		mod,
		master[[4]] %>%
			select(
				Season
			),
		permutations = 999)
	plot(ef)
}



{ # Overall Cross-Catchment Graphs
	
}

{ # DESeq2
	# Separate metadata and ASV count data
	asv_counts <- bact_collapsed[[4]] %>% select_after("invsimpson")
	meta <- bact_collapsed[[4]] %>% select(1:(which(names(bact_collapsed[[4]]) == "invsimpson")))
	
	# Round ASV counts (DESeq2 requires integer counts)
	asv_counts <- round(asv_counts)
	
	# Create DESeq2 dataset
	dds <- DESeqDataSetFromMatrix(
		countData = t(asv_counts),  # taxa in rows
		colData = meta,
		design = ~ Season  # change this to your grouping variable
	)
	
	# Filter out low-abundance ASVs
	dds <- dds[rowSums(counts(dds)) > 10, ]
	
	# Run DESeq2
	dds <- DESeq(dds, sfType = "poscounts")
	
	# Get results for group comparison (e.g. Summer vs Spring)
	res <- results(dds, contrast = c("Catchment", "Red", "Cober"))
	
	# View significant results
	res_sig <- res %>% as.data.frame() %>%
		rownames_to_column("Taxon") %>%
		filter(padj < 0.05)
	
	# Save results
	write_csv(res_sig, "deseq2_significant_taxa.csv")
}
