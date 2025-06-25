{ # User-definable parameters ----
	keyvals <- list(
		metals = c(
			"Aluminium",
			"Arsenic",
			"Barium",
			"Boron",
			"Cadmium",
			"Chromium",
			"Copper",
			"Iron",
			"Lead",
			"Lithium",
			"Magnesium",
			"Manganese",
			"Nickel",
			"Zinc"
		),
		# metals in mgL for conversion to μgL
		milligrams = c(
			"Sodium",
			"Potassium",
			"Magnesium",
			"Calcium"
		),
		subsampling_depth = 9884,
		catchments = c(
			"Hayle",
			"Cober",
			"Red",
			"Carnon"
		),
		a_indices = c(
			"richness",
			"shannon",
			"simpson",
			"invsimpson"
		),
		years = c(
			2013,
			2024
		),
	# number of taxonomic levels
		n_levels = 4,
	# number of samples sequenced
		# updated from 319 after new data received Jun 10th 2025
		n_samples_seq = 343,
	# maximum number of categories for numeric factors on MDS plots
		bins = 4,
	# the name of each taxonomic level in order
		# removed level 5; missing in the updated PCR data received on Jun 10th 2025
		PCR_levels = c(
			"Phylum",
			"Class",
			"Order",
			"Genus"
		)
	)
	# water quality parameters
		keyvals$parameters = paste0(
			keyvals$metals, ", Dissolved"
		)
	# 
		keyvals$medians = paste0(
			keyvals$metals, "_median.μgL"
		)
	# environmental factors of interest (including metals defined above)
		keyvals$com_factors = c(
			"Month"
		)
		keyvals$com_factors_proper = c(
			"Temp",
			"Month",
			"Total_Mean.mgL",
			keyvals$metals
		)
}
