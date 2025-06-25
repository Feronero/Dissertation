{
	{
		# master list of all paired EA points, bacteria points, and metal data
		master_paired <- list()
		# bact_points:metal
		for (i in 1:keyvals$n_levels) {
			master_paired[[i]] <- merge(
				x    = bact_data[ , c(1,2,3,4,5,6,15,16,17,18)],
				y    = metal_points$paired[ , c(4,10,14:ncol(metal_points$paired))],
				by.x = "Site",
				by.y = "eDNA_point"
			)
			# bact_points:metal:qPCR_data
			master$paired[[i]] <- merge(
				x = master$paired[[i]],
				y = qPCR_data$meaned,
				by = "Sample"
			)
			# bact_points:metal:qPCR_data:a-div
			master$paired[[i]] <- merge(
				x = master$paired[[i]],
				y = a_div[[i]],
				by = "SampleID"
			)
			# bact_points:metal:qPCR_data:a-div:PCR_data
			master$paired[[i]] <- merge(
				x = master$paired[[i]],
				y = PCR_data[[i]],
				by = "SampleID"
			)
		}
	}
	# only EA points with records of all metals ("complete records")
	{
		master_completemetals <- list()
		# bact_points:metal
		for (i in 1:keyvals$n_levels) {
			master$completemetals[[i]] <- merge(
				x    = bact_data$all[ , c(1,2,3,4,6,14,15,16,17,18)],
				y    = metal_points$complete[ , c(4,10,14:ncol(metal_points$complete))],
				by.x = "Site",
				by.y = "eDNA_point"
			)
			# bact_points:metal:qPCR_data
			master$completemetals[[i]] <- merge(
				x = master$completemetals[[i]],
				y = qPCR_data$meaned,
				by = "Sample"
			)
			# bact_points:metal:qPCR_data:a-div
			master$completemetals[[i]] <- merge(
				x = master$completemetals[[i]],
				y = a_div[[i]],
				by = "SampleID"
			)
			# bact_points:metal:qPCR_data:a-div:PCR_data
			master$completemetals[[i]] <- merge(
				x = master$completemetals[[i]],
				y = PCR_data[[i]],
				by = "SampleID"
			)
		}
	}
	# only EA points with complete records
	# only metals with variance
	{
		master$completewvariance <- list()
		# bact_points:metal
		for (i in 1:keyvals$n_levels) {
			master$completewvariance[[i]] <- merge(
				x = bact_data$all[ , c(1,2,3,4,6,14,15,16,17,18)],
				y = metal_points$complete_wvariance[ , c(4,10,14:ncol(metal_points$complete_wvariance))],
				by.x = "Site",
				by.y = "eDNA_point"
			)
			# bact_points:metal:qPCR_data
			master$completewvariance[[i]] <- merge(
				x = master$completewvariance[[i]],
				y = qPCR_data$meaned,
				by = "Sample"
			)
			# bact_points:metal:qPCR_data:a-div
			master$completewvariance[[i]] <- merge(
				x = master$completewvariance[[i]],
				y = a_div[[i]],
				by = "SampleID"
			)
			# bact_points:metal:qPCR_data:a-div:PCR_data
			master$completewvariance[[i]] <- merge(
				x = master$completewvariance[[i]],
				y = PCR_data[[i]],
				by = "SampleID"
			)
		}
	}
	# remove Arsenic from consideration to allow a few more sampling points to be used
	# remove sampling points until all catchments have 3 sampling points
	# some points have no sampling data from October, so remove October from consideration
	# Lithium now has zero variance, so needs to be excluded from consideration as well
	{
		master$balanced <- list()
		for (i in 1:keyvals$n_levels) {
			master$balanced[[i]] <- master$paired[[i]] %>%
				# begin filtering
				select(
					-Arsenic_median.μgL,
					-Boron_median.μgL,
					-Chromium_median.μgL,
					-Lithium_median.μgL
				) %>%
				filter(
					if_all(ends_with("median.μgL"), ~ !is.na(.x))
				) %>%
#				filter(!Month %in% c("October")) %>%
				# 1. Identify sites with data for all required months
				group_by(Catchment, Site) %>%
				filter(all(c("January", "April", "July", "October") %in% Month)) %>%
				ungroup()
			
			# 2. Randomly select 3 sites per river (from those with required monthly coverage)
			set.seed(123)  # For reproducibility
			selected_sites <- master$balanced[[i]] %>%
				group_by(Catchment) %>%
				distinct(Site) %>%
				slice_sample(n = 3) %>%   # no. sites per catchment
				ungroup()
			
			# 3. Retain 1 entry per month for selected sites
			master$balanced[[i]] <- master$balanced[[i]] %>%
				# Keep only selected sites
				semi_join(selected_sites, by = c("Catchment", "Site")) %>%
				group_by(Catchment, Site, Month) %>%
				slice_sample(n = 1) %>% # For each site-month, keep 1 random entry
				ungroup() %>%
				# Sort for clarity
				arrange(Catchment, Site, Month)
		}
	}
}

# Check for missing River -> Site -> Month combinations ----
{
	checklist <- list()
	checklist <- merge(
		x  = bact_data$all,
		y  = PCR_data[[1]],
		by = "SampleID"
	)
	
	# Check sample counts
	sample_counts <- checklist %>%
		group_by(Catchment, Site, Month) %>%
		summarise(Samples = n(), .groups = "drop")
	
	# View combinations with incorrect sample numbers
	incorrect_counts <- sample_counts %>% 
		filter(Samples != 3)
	
	# Print results
	if(nrow(incorrect_counts) == 0) {
		message("All site/month/river combinations have exactly 3 samples!")
	} else {
		message("Some combinations have incorrect sample counts:")
		print(incorrect_counts)
	}
	
	# Optional: Check total expected vs. observed samples
	expected_samples <- 29 * 4 * 3
	actual_samples <- nrow(checklist)
	cat("\nExpected samples:", expected_samples)
	cat("\nActual samples:", actual_samples)
	
	rm(checklist)
}
