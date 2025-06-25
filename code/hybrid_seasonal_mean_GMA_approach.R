## ── 0  Libraries ──────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(mgcv)      # GAMs
library(vegan)     # PERMANOVA / dbRDA

## ── 1  Read & prepare data ───────────────────────────────────────────────────
bact  <- bact_data[[4]]              # ← your already-rarefied OTU table
metal <- metal_wide                  # ← Environment Agency chemistry

## Bacterial metadata / OTU matrix (no rare-taxon filtering)
otu_cols <- grep("^d__", names(bact), value = TRUE)
meta     <- bact %>% mutate(Site = as.character(Site))
otu_rel  <- bact[otu_cols] %>% as.matrix() %>% sweep(1, rowSums(.), "/")
bcurt    <- vegdist(otu_rel, method = "bray")

## Metal table → long, with critical factor conversions
metal_long <- metal %>% 
	pivot_longer(Lead:Calcium, names_to = "metal",
					 values_to = "conc", values_drop_na = TRUE) %>% 
	mutate(
		season    = factor(Season,  levels = c("Winter","Spring","Summer","Autumn")),
		year      = factor(year),                 # RANDOM-EFFECT factor
		Site      = factor(Site),                 # RANDOM-EFFECT factor
		Catchment = factor(Catchment)
	)

## ── 2  Helper functions ──────────────────────────────────────────────────────
fit_gam <- function(df) {
	if (nrow(df) < 3 || var(df$conc) == 0) return(NULL)
	f <- if (n_distinct(df$season) > 1)
		conc ~ season + s(year, bs = "re") + s(Site, bs = "re")
	else  conc ~              s(year, bs = "re") + s(Site, bs = "re")
	tryCatch(gam(f, data = df, method = "REML"), error = function(e) NULL)
}

predict_gam <- function(mod) {
	if (is.null(mod)) return(NULL)
	expand_grid(
		Site   = levels(mod$model$Site),
		season = levels(mod$model$season),
		year   = levels(mod$model$year)
	) %>%
		mutate(fit = predict(mod, newdata = ., type = "response"))
}

## ── 3  Nest → model → predict ───────────────────────────────────────────────
metal_nested <- metal_long %>%
	group_by(Catchment, metal) %>%
	nest() %>%
	mutate(model = map(data, fit_gam),
			 pred  = map(model, predict_gam)) %>%
	select(-data, -model) %>%           # drop bulky columns
	unnest(pred)                        # long table of predictions

## ── 4  Attach predictions to bacterial metadata ─────────────────────────────
meta_enriched <- meta %>%
	mutate(year   = factor(year(Date)),            # match factor type
			 season = factor(Season,
			 					 levels = levels(metal_long$season))) %>%
	left_join(metal_nested %>%
				 	pivot_wider(names_from = metal,
				 					values_from = fit,
				 					names_glue = "{metal}_pred"),
				 by = c("Site","season","year","Catchment"))

## ── 5  Community response analysis (PERMANOVA shown) ────────────────────────
adonis2(
	bcurt ~ Lead_pred + Zinc_pred + Copper_pred + season,
	data   = meta_enriched,
	method = "bray",
	strata = meta_enriched$Site,   # replicates nested in Site
	by     = "margin"
)