#	install.packages("leaflet")
	library(leaflet)
	
	# Create a palette function based on catchment levels
	grouping_palette <- colorFactor(
		palette = hcl.colors(
			n = nlevels(filtered_results$catchment),  # Number of levels in catchment
			palette = "Dark 3"
		),
		domain = filtered_points$catchment  # Ensure this column exists in filtered_points
	)
	
	# Build the leaflet map
	map <- leaflet(
		data = filtered_points
#		options = leafletOptions(
#			minZoom = 10,
#			maxZoom = 100
#		)
	) %>%
		addTiles() %>%
		addCircleMarkers(
			# Dynamic marker coloring based on catchment
			lng = ~long,
			lat = ~lat,
			label = ~label,
			color = ~"black",
			fillColor = ~grouping_palette(catchment),
			fillOpacity = 1,
			stroke = TRUE
		) %>%
		addLegend(
			position = "bottomright",
			pal = grouping_palette,  # Use the palette function
			values = ~catchment,  # Reference the same column in the data
			title = "Catchment",
			labFormat = labelFormat(prefix = ""),
			opacity = 1
		)
	
	# Display the map
	map
