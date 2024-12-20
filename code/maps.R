#	install.packages("leaflet")
	library(leaflet)
	
	# Create a palette function based on catchment levels
	grouping_palette <- colorFactor(
		palette = hcl.colors(
			n = nlevels(filtered_points$catchment),  # Number of levels in catchment
			palette = "Dark 3"
		),
		domain = filtered_points$catchment  # Ensure this column exists in filtered_points
	)
	
	# Build the leaflet map
	map <- leaflet(
		data = filtered_points
	) %>%
		addTiles() %>%
		addCircleMarkers(
			# Dynamic marker colouring based on catchment
			lng = ~long,
			lat = ~lat,
			label = ~label,
			color = "black",
			fillColor = ~grouping_palette(catchment),
			fillOpacity = 1,
			stroke = TRUE
		) %>%
		addLegend(
			position = "bottomright",
			pal = grouping_palette,
			values = ~catchment,
			title = "Catchment",
			labFormat = labelFormat(prefix = ""),
			opacity = 1
		)
	map %>% addCircleMarkers(
		data = sample_metadata,
		lng = ~Longitude,
		lat = ~Latitude,
		label = ~Location,
		color = "black",
		fillColor = "white",
		fillOpacity = 0,
		stroke = TRUE
	)
	
	test_map <- leaflet(data = sample_metadata) %>%
		addTiles() %>%
		addCircleMarkers(
			lng = ~Longitude,
			lat = ~Latitude,
			label = ~Location,
			color = "black",
			fillColor = "white",
			fillOpacity = 1,
			stroke = TRUE
		)
	test_map

	map

	
	# Create a palette function based on catchment levels
	grouping_palette <- colorFactor(
		palette = hcl.colors(
			n = nlevels(filtered_points$catchment),  # Ensure levels are consistent
			palette = "Dark 3"
		),
		domain = filtered_points$catchment
	)
	
	# Build the leaflet map
	map <- leaflet(data = filtered_points) %>%
		addTiles() %>%
		addCircleMarkers(
			# Dynamic marker coloring based on catchment
			lng = ~long,
			lat = ~lat,
			label = ~label,
			color = "black",
			fillColor = ~grouping_palette(catchment),
			fillOpacity = 1,
			stroke = TRUE
		) %>%
		addLegend(
			position = "bottomright",
			pal = grouping_palette,
			values = ~catchment,
			title = "Catchment",
			labFormat = labelFormat(prefix = ""),
			opacity = 1
		) %>%
		# Add second layer of points
		addCircleMarkers(
			data = sample_metadata,
			lng = ~Longitude,
			lat = ~Latitude,
			label = ~Location,
			color = "black",
			fillColor = "white",
			fillOpacity = 1,
			stroke = TRUE
		)
	
	