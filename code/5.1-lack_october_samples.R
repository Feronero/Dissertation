### Many sites lack samples from October

set.seed(123)
# remove all bacteria and metal data from sites with incomplete metal records
balanced <- list()
balanced$all <- master[[5]] %>%
	filter(
		!notation %in% c(
			"SW-82221881", "SW-82310183",
			"SW-82310186", "SW-81950522", "SW-81950901"
		)
	) %>%
#				filter(
#					!Month %in% "October"
#				) %>%
	select(-Arsenic_mean.μgL) %>%
	group_by(Catchment, Site) %>%
	filter(all(c("January", "April", "July", "October") %in% Month)) %>%
	group_by(Catchment, Site, Month) %>%
	slice_sample(n = 1) %>%  # Or use slice_head(n = 1) for deterministic selection
	ungroup() %>%
	arrange(Catchment, Site, Month)  # Optional: Sort for clarity

View(balanced$all)
