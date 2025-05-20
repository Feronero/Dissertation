library(geosphere)  # For distGeo()

# Example data: Replace with actual dataset
bacteria_points <- data.frame(
	id = 1:5,
	lon = c(-4.5, -4.6, -4.7, -4.55, -4.65),
	lat = c(50.2, 50.25, 50.3, 50.22, 50.27)
)

metal_points <- data.frame(
	id = 101:105,
	lon = c(-4.51, -4.61, -4.71, -4.52, -4.63),
	lat = c(50.21, 50.26, 50.31, 50.23, 50.28),
	order = c(1, 2, 3, 4, 5)  # Order from source to mouth
)

# Define function to compute geospatial distance
compute_distance <- function(lat1, lon1, lat2, lon2) {
	return(distGeo(c(lon1, lat1), c(lon2, lat2)))  # Returns meters
}

# Define function to determine upstream/downstream based on metal order
is_upstream <- function(bacteria_point, metal_point) {
	return(metal_point$order < bacteria_point$nearest_metal_order)  # Lower order = upstream
}

# Define Influence Index function
calculate_influence_index <- function(distance, upstream, k = 0.2) {
	W_u <- ifelse(upstream, 1, k)  # Upstream = 1, Downstream = k (penalty)
	return(W_u / (distance + 1))
}

# Maximum pairing radius (in meters)
max_distance <- 500

# Find nearest metal sampling point (by distance) for each bacteria point
bacteria_points$nearest_metal_order <- NA

for (i in 1:nrow(bacteria_points)) {
	b <- bacteria_points[i, ]
	distances <- sapply(1:nrow(metal_points), function(j) {
		compute_distance(b$lat, b$lon, metal_points$lat[j], metal_points$lon[j])
	})
	in_range <- which(distances <= max_distance)
	if (length(in_range) > 0) {
		nearest_metal <- metal_points[in_range[which.min(distances[in_range])], ]
		bacteria_points$nearest_metal_order[i] <- nearest_metal$order
	}
}

# Pairing process
pairings <- data.frame(bacteria_id = integer(), metal_id = integer(), distance = numeric(), influence_index = numeric())

for (i in 1:nrow(bacteria_points)) {
	b <- bacteria_points[i, ]
	if (is.na(b$nearest_metal_order)) next  # Skip if no metal point within range
	
	candidates <- data.frame()
	
	for (j in 1:nrow(metal_points)) {
		m <- metal_points[j, ]
		dist <- compute_distance(b$lat, b$lon, m$lat, m$lon)
		
		if (dist <= max_distance) {
			upstream <- is_upstream(b, m)
			II <- calculate_influence_index(dist, upstream)
			candidates <- rbind(candidates, data.frame(metal_id = m$id, distance = dist, influence_index = II))
		}
	}
	
	if (nrow(candidates) > 0) {
		best_match <- candidates[which.max(candidates$influence_index), ]
		pairings <- rbind(pairings, data.frame(
			bacteria_id = b$id,
			metal_id = best_match$metal_id,
			distance = best_match$distance,
			influence_index = best_match$influence_index
		))
	}
}

# View the final pairings
print(pairings)
