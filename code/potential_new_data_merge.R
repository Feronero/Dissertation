# Master data compilation function
compile_master_data <- function(bact_data, metal_data, qPCR_data, a_div, PCR_data, 
										  level, metal_cols = NULL, filter_complete = FALSE, 
										  filter_variance = FALSE) {
	# Start with basic merge of bacterial and metal data
	df <- merge(
		x = bact_data[, c(1, 2, 3, 4, 6, 15, 16, 17, 18, 20)],
		y = metal_data[, c(4, 10, 14:ncol(metal_data))],
		by.x = "Site",
		by.y = "eDNA_point"
	)
	
	# Apply metal filters if specified
	if (!is.null(metal_cols)) {
		df <- df %>% select(all_of(metal_cols))
	}
	if (filter_complete) {
		df <- df %>% filter(if_all(ends_with("median.μgL"), ~!is.na(.x)))
	}
	
	# Merge with remaining datasets
	df <- df %>%
		merge(y = qPCR_data$meaned, by = "Sample") %>%
		merge(y = a_div[[level]], by = "SampleID") %>%
		merge(y = PCR_data[[level]], by = "SampleID")
	
	return(df)
}

# Balanced sampling function
create_balanced_subsample <- function(df, sites_per_catchment = 3, months = c("January", "April", "July", "October")) {
	# Filter to sites with data for all required months
	complete_sites <- df %>%
		group_by(Catchment, Site) %>%
		filter(all(months %in% Month)) %>%
		ungroup()
	
	# Randomly select sites per catchment
	set.seed(123)  # For reproducibility
	selected_sites <- complete_sites %>%
		group_by(Catchment) %>%
		distinct(Site) %>%
		slice_sample(n = sites_per_catchment) %>%
		ungroup()
	
	# Select one sample per site-month combination
	balanced_data <- complete_sites %>%
		semi_join(selected_sites, by = c("Catchment", "Site")) %>%
		group_by(Catchment, Site, Month) %>%
		slice_sample(n = 1) %>%
		ungroup() %>%
		arrange(Catchment, Site, Month)
	
	return(balanced_data)
}

# Main data compilation
master <- list(
	# All paired points with minimal filtering
	paired = lapply(1:keyvals$n_levels, function(i) {
		compile_master_data(bact_data, metal_points$paired, qPCR_data, a_div, PCR_data, i)
	}),
	
	# Only points with complete metal records
	completemetals = lapply(1:keyvals$n_levels, function(i) {
		compile_master_data(bact_data, metal_points$complete, qPCR_data, a_div, PCR_data, i, 
								  filter_complete = TRUE)
	}),
	
	# Only metals with variance
	completewvariance = lapply(1:keyvals$n_levels, function(i) {
		compile_master_data(bact_data, metal_points$complete_wvariance, qPCR_data, a_div, PCR_data, i,
								  filter_complete = TRUE)
	}),
	
	# Balanced design subset
	balanced = lapply(1:keyvals$n_levels, function(i) {
		compile_master_data(bact_data, metal_points$paired, qPCR_data, a_div, PCR_data, i,
								  metal_cols = c(-Arsenic_median.μgL, -Boron_median.μgL, 
								  					-Chromium_median.μgL, -Lithium_median.μgL),
								  filter_complete = TRUE) %>%
			create_balanced_subsample()
	})
)

# Data quality check function
check_sample_balance <- function(data) {
	sample_counts <- data %>%
		group_by(Catchment, Site, Month) %>%
		summarise(Samples = n(), .groups = "drop")
	
	incorrect_counts <- sample_counts %>% 
		filter(Samples != 3)
	
	if (nrow(incorrect_counts) == 0) {
		message("All site/month/river combinations have exactly 3 samples!")
	} else {
		message("Some combinations have incorrect sample counts:")
		print(incorrect_counts)
	}
	
	expected_samples <- 29 * 4 * 3
	actual_samples <- nrow(data)
	cat("\nExpected samples:", expected_samples)
	cat("\nActual samples:", actual_samples)
}

# Run quality check on first level data
check_sample_balance(
	merge(
		x = bact_data$all,
		y = PCR_data[[1]],
		by = "SampleID"
	)
)