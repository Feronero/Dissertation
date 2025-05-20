# bar chart ----
{
	bact_col_names <- 
	
	df_agg <- master$completewvariance[[1]] %>%
		group_by(Catchment, Month) %>%
		summarise(across(all_of(colnames)))
		summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
		ungroup()
	top_phyla <- df_agg %>%
		summarise(across(-c(Catchment, Month), mean)) %>%
		pivot_longer(everything(), names_to = "phylum", values_to = "abundance") %>%
		arrange(desc(abundance)) %>%
		slice_head(n = 10) %>%
		pull(phylum)
	
	# Combine non-top phyla into "Other"
	df_agg_top <- df_agg %>%
		mutate(Other = rowSums(across(!all_of(c("Catchment", "Month", top_phyla))))) %>%
		dplyr::select(all_of(c("Catchment", "Month", top_phyla, "Other"))) %>%
		pivot_longer(-c(Catchment, Month), names_to = "phylum", values_to = "abundance")
	# Order seasons chronologically
	season_order <- c("January", "April", "July", "October")
	df_agg_top$Month <- factor(df_agg_top$Month, levels = season_order)
	
	# Create the plot
	ggplot(df_agg_top, aes(x = interaction(Catchment, Month, sep = " - "), y = abundance, fill = phylum)) +
		geom_col(position = "stack", width = 0.8) +
		scale_fill_viridis_d(option = "D") + # Use a colorblind-friendly palette
		scale_x_discrete(
			breaks = c("River1 - January", "River2 - January", "River3 - January", "River4 - January"),
			labels = c("River 1", "River 2", "River 3", "River 4")
		) +
		labs(
			x = "River (grouped by season: Jan, Apr, Jul, Oct)",
			y = "Relative Abundance (%)",
			fill = "Phylum",
			title = "Top 10 Bacterial Phyla by River and Season"
		) +
		theme_minimal() +
		theme(
			axis.text.x = element_text(angle = 45, hjust = 1),
			panel.spacing = unit(1, "cm"), # Add space between river clusters
			legend.position = "bottom"
		)
}
