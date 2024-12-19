	install.packages("leaflet")
	library(leaflet)

map <- leaflet(
	data = ,
	options = leafletOptions(
		minZoom = 10,
		maxZoom = 100)
) |> 
	
	#fitBounds(
	
	#) %>%
	
#	setView(
#		lng =
#			min(patches$long) +
#		lat =
#			(diff(range(patches$long))/2),
#			min(patches$lat) +
#			(diff(range(patches$lat))/2),
#		zoom = 10
#	) %>% 
	
	addTiles(
	) %>%
	
	addCircleMarkers(
		#clusterOptions = markerClusterOptions(),
		lng = ~long,
		lat = ~lat,
		label = ~number,
		color = ~"black",
		fillColor = "white",
		fillOpacity = 1,
		stroke = TRUE
	) %>%
	
	addLegend(
		position = "bottomright",
		pal = shades_degree,
		values = ~degree,
		title = "Degree",
		labFormat = labelFormat(prefix = ""),
		opacity = 1
	)

for (poly in polylines){
	map <- map %>%
		addPolylines(
			lat = poly$lat,
			lng = poly$long,
			color = "black"
		)
}

map