desc <- list(
	a_div = "richness, shannon, simpson, and inverse simpson of all samples separated by taxonomic level",
	bact_data = "environmental metadata of all bacteria samples",
	bact_points = "list of all bacteria sampling points",
	keyvals = "user-defined terms specific for this investigation",
	metal_data = list(
		all = "all water quality parameters for each of the catchments from all 10 years",
		paired = "all water quality parameters for all paired metal points",
		filtered = "metal parameters of interest for all paired metal points",
		reference = "the units of each metal parameter of interest",
		wide = "wide-format metal data of interest; each metal has its own column"
	),
	metal_points = list(
		all = "location metadata of all metal sampling points. No water quality parameters included",
		paired = "location and metal data of interest of all paired metal points",
		complete = "location + metal data from metal points with >1 measurement of each metal",
		complete_wvariance = "location + metal data from metal points with >1 measurement of each metal WITHOUT metal that doesn't vary (some metals seem to be held at a constant level across samples and rivers)"
	),
	PCR_data = "sequence count data, seperated by taxonomic level",
	qPCR_data = "relative starting overall bacterial abundance. Not absolute; no standard curve used"
)
