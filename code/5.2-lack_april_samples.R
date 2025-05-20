### Samples lack either Oct, Jan, or July

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
					filter(
						!Month %in% c("October", "January")
					) %>%
	select(-Arsenic_mean.μgL) %>%
	group_by(Catchment, Site) %>%
	filter(all(c("July", "April") %in% Month)) %>%
	group_by(Catchment, Site, Month) %>%
	slice_sample(n = 1) %>%  # Or use slice_head(n = 1) for deterministic selection
	ungroup() %>%
	arrange(Catchment, Site, Month)

View(balanced$all)

###### The best I can manage - if any varibles are increased, becomes unbalanced

balanced <- list()
balanced$all <- master[[5]] %>%
	filter(
		!notation %in% c(
			"SW-82221881", "SW-81950901",
			"SW-82310186", "SW-81950522"
		)
	) %>%
	filter(
		!Month %in% c("October")
	) %>%
	select(-Arsenic_mean.μgL) %>%
# 1. Identify sites with data for all 4 months (Jan-Apr)
	group_by(Catchment, Site) %>%
	filter(all(c("January", "April", "July") %in% Month)) %>%
	ungroup()

# 2. Randomly select 3 sites per river (from those with full monthly coverage)
set.seed(123)  # For reproducibility
selected_sites <- balanced$all %>%
	group_by(Catchment) %>%
	distinct(Site) %>%  # Unique sites per river
	slice_sample(n = 3) %>%    # Cap at 3 sites; selects fewer if insufficient
	ungroup()

# 3. Retain 1 entry per month for selected sites
balanced$all <- balanced$all %>%
	# Keep only selected sites
	semi_join(selected_sites, by = c("Catchment", "Site")) %>%
	# For each site-month, keep 1 entry (randomly selected)
	group_by(Catchment, Site, Month) %>%
	slice_sample(n = 1) %>%
	ungroup() %>%
	# Sort for clarity
	arrange(Catchment, Site, Month)

View(balanced$all)

