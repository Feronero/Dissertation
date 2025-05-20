par(
	mfrow = c(1,1)
)
hist(
	data$richness
)
	par(
		mfrow = c(2,2)
	)
	plot(
		richness ~ PC1,
		data = data,
		main = "PC1"
	)
	plot(
		richness ~ Temp,
		data = data,
		main = "Temp"
	)
	plot(
		richness ~ pH,
		data = data,
		main = "pH"
	)
	plot(
		richness ~ PC2,
		data = data,
		main = "PC2"
	)
	par(
		mfrow = c(1,1)
	)
	lm1 <- glm(
		formula = richness ~ PC1 + PC2 + PC3 + PC4 + Temp + Month + pH + DTM.km,
		family = "gaussian" (link = "identity"),
		data = data
	)
	plm1 <- glm(
		formula = richness ~ PC1 + PC2 + PC3 + PC4 + Temp + Month + pH,
		family = "poisson" (link = "log"),
		data = data
	)
	bc <- boxcox(plm1)
	plot(lm1)
	
	lambda <- bc$x[which.max(bc$y)]
	plm2 <- glm(
		formula = ((richness^lambda-1)/lambda) ~ PC1 + PC2 + PC3 + PC4 + Temp + Month + pH,
		family = "poisson" (link = "log"),
		data = data
	)
	nb_mod <- glm.nb(richness ~ Temp + pH + Altitude + Month + PC1 + PC2 + PC3 + PC4, data = data)
	summary(nb_mod)
	gam_mod <- gam(richness ~ s(Temp) + s(pH) + Month + s(PC1) + PC2 + PC3 + PC4, 
				family = poisson(link = "log"), data = data)
	
	summary(gam_mod)
	plot(gam_mod, page = 1)
	
	
	
	plm3 <- glm(
		formula = richness ~ Temp + PC1 + Month * Catchment,
		family = "poisson" (link = "log"),
		data = data
	)
	par(mfrow=c(2,2))
	plot(plm3)
	par(mfrow=c(1,1))
	
#	Check for cross-catchment effect
#	
	