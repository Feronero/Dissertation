query <- unique(
	bact_points$all$Site
)
grepl(
	pattern = paste0(
		query,
		collapse = "|"
	),
)