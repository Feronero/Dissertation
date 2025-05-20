#	install.packages("leaflet")
	library(leaflet)

# Create Map
{
	# temporary points storage
	unique <- bact_points$all$Site
	# Create a palette function based on catchment levels
	grouping_palette <- colorFactor(
		palette = hcl.colors(
			n = nlevels(metal_points$all$catchment),  # Ensure levels are consistent
			palette = "Dark 3"
		),
		domain = metal_points$all$catchment
	)
	
	# Build the leaflet map
	map <- leaflet(data = metal_points$all) %>%
		addTiles() %>%
		addCircleMarkers(
			# Dynamic marker coloring based on catchment
			lng = ~long,
			lat = ~lat,
			label = ~Name,
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
			data = unique,
			lng = ~Longitude,
			lat = ~Latitude,
			label = ~Site,
			color = "black",
			fillColor = "white",
			fillOpacity = 0.2,
			stroke = TRUE
		)
	map
}
	