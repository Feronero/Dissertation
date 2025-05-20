{
	b_div <- list()
	for (j in keyvals$catchments) {
		for (x in keyvals$com_factors) {
			for (i in 1:keyvals$n_levels) {
				# PCR and sampling point data for all catchments across all tax. levels
				b_div[[j]]$data[[i]] <- subset(
					x = bact_master[[i]],
					subset = Catchment == j
				)
				
				# rename community data selection parameters for convenience
				community_data <- b_div[[j]]$data[[i]][, 24:ncol(b_div[[j]]$data[[i]])]
				community_factors <- as.factor(b_div[[j]]$data[[i]][, x])
				
				# MDS results for all tax. levels
				b_div[[j]]$MDS[[i]] <- metaMDS(
					comm = community_data,
					distance = "bray",
					k = 2,
					trymax = 100
				)
				# ADONIS results for all factors and tax.levels
				b_div[[j]]$ADONIS[[x]][[i]] <- adonis2(
					formula = community_data ~ community_factors
				)
				# plot and hulls for all factors and tax. levels
				b_div[[j]]$plot[[x]][[i]] <- plot(
					b_div[[j]]$MDS[[i]],
					type = "n",
					main = paste(
						j, x, i
					)
				)
				points(
					b_div[[j]]$MDS[[i]],
					pch = 16,
					col = as.numeric(
						community_factors	
					)
				)
				for (z in 1:nlevels(community_factors)) {
					# Subset points for each group
					group_points <- b_div[[j]]$MDS[[i]]$points[community_factors == levels(community_factors)[z], ]
					
					# Calculate convex hull
					hull <- chull(group_points)
					
					# Add polygon to plot
					polygon(group_points[hull, ], 
						   border = z,
						   lty = 2)
				}
			}
		}
	}
	rm(community_data)
}