# Set up data ----
{
	# dataframe of all metal points with complete data
	metal_points <- metal_points %>%
		filter(if_all(14:30, ~ !is.na(.x)))
	na_metals <- norm_tests$medianmetal %>%
		filter(is.na(W)) %>%
		pull(Metal)
	metal_points$complete_wvariance <- metal_points$complete %>%
		select(-any_of(na_metals))
	varying_metal_cols <- metal_points$complete_wvariance[, names(metal_points$complete_wvariance) %in% keyvals$medians]
	rm(na_metals)
}

# Non-Parametric Correlation Test on Metals ----
{
	# run corr test
	metal_data$cor_matrix <- as.data.frame(
		cor(
			varying_metal_cols,
			use    = "pairwise.complete.obs",
			method = "spearman"
		)
	)
	# export correlation table
	write_xlsx(
		x    = metal_data$cor_matrix,
		path = "outputs/tables/metal_cor.xlsx"
	)
	rm(varying_metal_cols)
}
