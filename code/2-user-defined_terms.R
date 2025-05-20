# User-definable parameters
{
	keyvals <- list(
		metals = c(
			"Aluminium",
			"Arsenic",
			"Barium",
			"Boron",
			"Cadmium",
			"Calcium",
			"Chromium",
			"Copper",
			"Iron",
			"Lead",
			"Lithium",
			"Magnesium",
			"Manganese",
			"Nickel",
			'Potassium',
			"Sodium",
			"Zinc"
		),
		# metals in mgL for conversion to μgL
		milligrams = c(
			"Sodium",
			"Potassium",
			"Magnesium",
			"Calcium"
		),
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
		n_levels = 5,
		# number of samples sequenced
		n_samples_seq = 319,
		# maximum number of categories for numeric factors on MDS plots
		bins = 4,
		# the name of each taxonomic level in order
		PCR_levels = c(
			"level_2",
			"level_3",
			"level_4",
			"level_5",
			"level_6"
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
