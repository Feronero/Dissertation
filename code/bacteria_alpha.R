{ # Trial models for testing PC1 ----
	# valid
		random2 <- lmer(
			PC1 ~
				pH + Temp + River_Gradient + Season +
				(1 | Catchment / Site),
			data = test[[4]]
		)
		check_model(random2)
		simulateResiduals(random2, plot = TRUE)
	# !! singular
		random2_int <- lmer(
			PC1 ~
				(pH + Temp + River_Gradient) * Season +
				(1 | Catchment / Site),
			data = test[[4]]
		)
		check_model(random2)
		simulateResiduals(random2, plot = TRUE)
	# valid
		random1 <- lmer(
			PC1 ~
				pH + Temp + River_Gradient + Season +
				(1 | Catchment),
			data = test[[4]]
		)
		check_model(random1)
		simulateResiduals(random1, plot = TRUE)
	# valid
		random1_int <- lmer(
			PC1 ~
				(pH + Temp + River_Gradient) * Season +
				(1 | Catchment),
			data = test[[4]]
		)
		check_model(random1_int)
		simulateResiduals(random1_int, plot = TRUE)
	# !! singular
		random1_int_3 <- lmer(
			PC1 ~ 
				(pH + Temp + River_Gradient) * Season +
			   (1 + Temp + pH + River_Gradient || Catchment),
			data = test[[4]]
		)
		VarCorr(random1_int_3)
	# MAYBE OKAY????? SKSDVOSDGVKFGLK;FGKDBFAKL;ABLKSGNK;;OSKNKO;K;DRG;KGRNGRH;KGLNR
		random1_int2 <- lmer(
			PC1 ~ 
				(pH + Temp + River_Gradient) * Season +
			   (1 + Temp + pH || Catchment),
			data = test[[4]]
		)
		check_model(random1_int2)
		VarCorr(random1_int2)
		simulateResiduals(random1_int2, plot = TRUE)
	# valid
		random1_int1 <- lmer(
			PC1 ~ 
				(pH + Temp + River_Gradient) * Season +
			   (1 + Temp || Catchment),
			data = test[[4]]
		)
		VarCorr(random1_int1)
		check_model(random1_int1)
		simulateResiduals(random1_int1, plot = TRUE)
	# valid
		fixed <- lm(
			PC1 ~
				pH + Temp + River_Gradient + Season + Catchment,
			data = test[[4]]
		)
		check_model(fixed)
		simulateResiduals(fixed, plot = TRUE)
	# choose model
		compare_performance(random2, random1, random1_int, random1_int1, fixed, rank = TRUE)
		test_performance(random2, random1, random1_int, random1_int1, fixed)
		simulateResiduals(toggle_random1_int, plot = TRUE)
		testOutliers(toggle_random1_int, plot = TRUE)
		report(random1)
}

{ # Trial models for testing PC2 ----
	# singular
		random2 <- lmer(
			PC2 ~
				pH + Temp + River_Gradient + Season +
				(1 | Catchment / Site),
			data = test[[4]]
		)
	# singular
		random1 <- lmer(
			PC2 ~
				pH + Temp + River_Gradient + Season +
				(1 | Catchment),
			data = test[[4]]
		)
	# 
		fixed <- lm(
			PC2 ~
				pH + Temp + River_Gradient + Season + Catchment,
			data = test[[4]]
		)
		check_model(fixed)
	#
		fixed_int <- lm(
			PC2 ~
				(pH + Temp + River_Gradient + Season) * Catchment,
			data = test[[4]]
		)
		check_model(fixed_int)
	# compare models	
		compare_performance(fixed, fixed_int)
		test_performance(fixed, fixed_int)
	# results
		report(fixed)
}

{ # Trial models for testing Cq_Median ----
	# valid
		random2 <- lmer(
			log(Cq_Median) ~
				pH + Temp + River_Gradient + Season +
				(1 | Catchment / Site),
			data = test[[4]]
		)
		check_model(random2)
		simulateResiduals(random2, plot = TRUE)
		random2 %>%
			simulateResiduals(n = 1000) %>%
			testOutliers(type = "bootstrap", nBoot = 1000) # solved
	# valid
		random2_int <- lmer(
			log(Cq_Median) ~
				(pH + Temp + River_Gradient) * Season +
				(1 | Catchment / Site),
			data = test[[4]]
		)
		check_model(random2_int)
		simulateResiduals(random2_int, plot = TRUE)
		random2_int %>%
			simulateResiduals(n = 1000) %>%
			testOutliers(type = "bootstrap", nBoot = 1000) # solved
	# valid
		random1 <- lmer(
			log(Cq_Median) ~
				pH + Temp + River_Gradient + Season +
				(1 | Catchment),
			data = test[[4]]
		)
		check_model(random1)
		simulateResiduals(random1, plot = TRUE) # sig outliers
		random1 %>%
			simulateResiduals(n = 1000) %>%
			testOutliers(type = "bootstrap", nBoot = 1000) # solved
	# valid
		random1_int <- lmer(
			log(Cq_Median) ~
				(pH + Temp + River_Gradient) * Season +
				(1 | Catchment),
			data = test[[4]]
		)
		check_model(random1_int)
		simulateResiduals(random1_int, plot = TRUE) # no sig outliers
	# !! singular - pH singular
		random1_3 <- lmer(
			log(Cq_Median) ~
				(pH + Temp + River_Gradient) * Season +
				(1 + Temp + River_Gradient + pH || Catchment),
			data = test[[4]]
		)
		check_model(random1_int)
		simulateResiduals(random1_int, plot = TRUE) # 
		VarCorr(random1_3)
	# valid
		random1_2 <- lmer(
			log(Cq_Median) ~
				(pH + Temp + River_Gradient) * Season +
				(1 + Temp + River_Gradient || Catchment),
			data = test[[4]]
		)
		check_model(random1_2)
		simulateResiduals(random1_int, plot = TRUE) # 
		VarCorr(random1_2)		
	# Choose model
		compare_performance(random2, random2_int, random1, random1_int, random1_2)
		test_performance(random2, random2_int, random1, random1_int, random1_2)
	# random2 best overall
		summary(random2)
		report(random2)
}

{ # Plot PC1 effect ----
	# calculate mean values
means <- test[[4]] %>% summarize(
  pH = mean(pH),
  Temp = mean(Temp),
  River_Gradient = mean(River_Gradient)
)

# create new data frame for Season levels
newdat <- expand.grid(
  pH = means$pH,
  Temp = means$Temp,
  River_Gradient = means$River_Gradient,
  Season = factor(c("Spring","Summer","Autumn","Winter"), levels = levels(test[[4]]$Season))
)

pred <- predict(random1, newdata = newdat,
                re.form = NA,  # no random intercepts
                se.fit = TRUE,
                type = "response") # predicting PC1

newdat$fit <- pred$fit
newdat$se <- pred$se.fit

newdat <- newdat %>%
  mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

ggplot(newdat, aes(x = Season, y = fit)) +

  # Layer 1: all seasons in grey
  geom_point(data = newdat,
             color = "grey50", size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat,
                aes(ymin = lower, ymax = upper),
                color = "grey50", width = 0.2, alpha = 1, show.legend = FALSE) +

  # Layer 2: overlay Winter in black
  geom_point(data = subset(newdat, Season == "Winter"),
             color = "black", size = 3, alpha = 1) +
  geom_errorbar(data = subset(newdat, Season == "Winter"),
                aes(ymin = lower, ymax = upper),
                color = "black", width = 0.2, alpha = 1) +

  # labs + theme
  labs(
    y = "Predicted PC1 Value",
    title = "Predicted PC1 by Season (holding covariates at mean)"
  ) +
  theme_classic()
}

{ # Plot Cq_Median effect ----
# calculate mean values
means <- test[[4]] %>% summarize(
  pH = mean(pH),
  Temp = mean(Temp),
  River_Gradient = mean(River_Gradient)
)

# create new data frame for Season levels
newdat <- expand.grid(
  pH = means$pH,
  Temp = means$Temp,
  River_Gradient = means$River_Gradient,
  Season = factor(c("Spring","Summer","Autumn","Winter"), levels = levels(test[[4]]$Season))
)

pred <- predict(random2, newdata = newdat,
                re.form = NA,  # no random intercepts
                se.fit = TRUE,
                type = "link") # predicting log(Cq_Median)

newdat$fit <- pred$fit
newdat$se <- pred$se.fit

newdat <- newdat %>%
  mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

newdat <- newdat %>%
  mutate(
    fit_resp   = exp(fit),
    lower_resp = exp(lower),
    upper_resp = exp(upper)
  )

ggplot(newdat, aes(x = Season, y = fit_resp)) +

  # Layer 1: all seasons in grey
  geom_point(data = newdat,
             color = "grey50", size = 3, alpha = 1, show.legend = FALSE) +
  geom_errorbar(data = newdat,
                aes(ymin = lower_resp, ymax = upper_resp),
                color = "grey50", width = 0.2, alpha = 1, show.legend = FALSE) +

  # Layer 2: overlay Autumn in black
  geom_point(data = subset(newdat, Season == "Autumn"),
             color = "black", size = 3, alpha = 1) +
  geom_errorbar(data = subset(newdat, Season == "Autumn"),
                aes(ymin = lower_resp, ymax = upper_resp),
                color = "black", width = 0.2, alpha = 1) +

  # labs + theme
  labs(
    y = "Predicted Median Cq Value",
    title = "Predicted log(Cq_Median) by Season (holding covariates at mean)"
  ) +
  theme_classic()
}

{ # Test plot - not valid without checks ----
# Convert data to long format so each predictor is in its own panel
d_long <- test[[1]] %>%
  pivot_longer(
    cols = c(pH, Temp, River_Gradient),
    names_to = "Predictor",
    values_to = "Value"
  )

# Plot with facets for each predictor and season-specific regression lines
ggplot(d_long, aes(x = Value, y = PC2, color = Season)) +
  facet_wrap(~ Predictor, scales = "free_x") +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Predictor value",
    y = "PC2",
    color = "Season",
    title = "Effect of pH, Temp, and River Gradient on PC2 by Season"
  )

check_model(toggle)
}
