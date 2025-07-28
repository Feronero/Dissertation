# 3. Create mvabund object
Yv <- mvabund(
	bact_collapsed[[4]] %>%
		drop_na(Cq_Median) %>%
		select(richness, Cq_Median, invsimpson) %>%
		scale()
)

# 4. Fit multivariate GLM (choose family; Gaussian is OK for continuous diversity)
mod <- manylm(
	Yv ~
		Temp + pH + River_Gradient,
   data = bact_collapsed[[4]] %>% drop_na(Cq_Median))

# 5. Check diagnostics: residual vs. fitted plots
plot(mod)  # looks for non-linearity, dispersion issues  [oai_citation:0‡RDocumentation](https://www.rdocumentation.org/packages/mvabund/versions/3.6.3/topics/manyglm) [oai_citation:1‡cran.r-universe.dev](https://cran.r-universe.dev/mvabund/doc/manual.html)

mod_gamma <- manylm(Yv ~ Temp + pH + River_Gradient + Season,
                     data = bact_collapsed[[4]] %>% drop_na(Cq_Median), family = "gamma")
plot.manylm(mod_gamma)

# 6. Test effect of predictors on community and each metric
anova_res <- anova(mod,
                   resamp = "pit.trap",
                   test = "LR",
                   p.uni = "unadjusted",
                   nBoot = 999)

print(anova_res)  # displays community-level and per-metric p-values  [oai_citation:2‡RDocumentation](https://www.rdocumentation.org/packages/mvabund/versions/4.2.1/topics/anova.manyglm)

# 7. Inspect coefficients for directionality
summary(mod)   # shows coefficients per metric (positive = increasing)
coefplot(mod)  # visualize across diversity metrics