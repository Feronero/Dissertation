nrow(subset(sample_metadata, Catchment == "Hayle"))

test <- subset(sample_metadata, Catchment == "Hayle")

length(unique(test$Site))

length(unique(subset(sample_metadata, Catchment == "Red")$Site))

query <- matrix(0, 29, 2)
for (i in 1:29) {
	float <- subset(
		x = sample_metadata,
		subset = Site == i
	)
	query[i, 1] <- nrow(
		float
	)
	query[i, 2] <- length(
		unique(
			float$Month
		)
	)
}
unique(sample_metadata$Month)
