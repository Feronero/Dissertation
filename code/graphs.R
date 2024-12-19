boxplot(
	formula = result.Copper_ug.l ~ notation,
	data = subset(
		x = filtered_results,
		subset = catchment == "hayle"
	)
)
