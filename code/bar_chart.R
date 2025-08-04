## bar chart ----
{
	## select bacterial columns
	bact_cols <- bact_data[[1]] %>%
		select_after("invsimpson") %>%
		names()
	bact_cols <- grep(
		pattern = "^d__",
		x = names(bact_data[[1]]),
		value = TRUE
	)
	## converts raw sequence counts to relative abundance percentages
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
		dplyr::select(
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
	season_order <- c("Spring", "Summer", "Autumn", "Winter")
	df_agg_top$Season <- factor(
		x = df_agg_top$Season,
		levels = season_order
	)
	
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
