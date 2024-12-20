for (i in 1:length(catchments)) {
	boxplot(
		formula = result.Copper_ug.l ~ notation,
		main = paste0(catchments[i]),
		data = subset(
			x = filtered_results,
			subset = catchment == catchments[i]
		)
	)
}
