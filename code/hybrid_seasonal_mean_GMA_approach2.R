## ──────────────────────────────────────────────────────────────────────────────
##  Hybrid metal-exposure pipeline
##  Option A: seasonal mean (all sites)
##  Option B: year-smooth GAM  (high-resolution sites only, ≥ 40 records)
## ──────────────────────────────────────────────────────────────────────────────
library(tidyverse)
library(lubridate)
library(mgcv)
library(vegan)          # for Bray–Curtis + adonis2

## ── 1  Read data ─────────────────────────────────────────────────────────────
metal <- metal_wide
bact  <- bact_data[[4]]     # already rarefied

## ── 2  Bacterial metadata & distance matrix (unchanged) ──────────────────────
otu_cols <- grep("^d__", names(bact), value = TRUE)
meta <- bact %>%
	mutate(Site      = as.character(Site),
			 Catchment = factor(Catchment),
			 season    = factor(Season, levels = c("Winter","Spring","Summer","Autumn")),
			 year      = factor(year(Date)))          # for later joins

otu_rel <- bact[otu_cols] %>% as.matrix() %>% sweep(1, rowSums(.), "/")
bcurt   <- vegdist(otu_rel, method = "bray")

## ── 3  Metal table → long & tidy ─────────────────────────────────────────────
metal_long <- metal %>%
	pivot_longer(Lead:Calcium, names_to = "metal", values_to = "conc",
					 values_drop_na = TRUE) %>%
	mutate(season    = factor(Season,  levels = levels(meta$season)),
			 year      = factor(year),
			 Site      = factor(Site),
			 Catchment = factor(Catchment)) %>%
	select(Site, Catchment, year, season, metal, conc)

## ── 4  Identify high-resolution sites (≥ 40 rows) ───────────────────────────
hi_res_sites <- metal_long %>%
	count(Site, name = "n") %>%
	filter(n >= 40) %>%
	pull(Site)

## ── 5  Option A  (seasonal means, all sites)──────────────────────────────────
pred_A <- metal_long %>%
	group_by(Site, Catchment, season, metal) %>%
	summarise(fit = mean(conc, na.rm = TRUE), .groups = "drop") %>%
	pivot_wider(names_from = metal, values_from = fit,
					names_glue = "{metal}_pred")

## ── 6  Option B  (GAM for hi-res sites)───────────────────────────────────────
campaign_years <- meta$year |> levels()            # bacterial sampling years

safe_gam <- function(df) {                         # numeric smooth on year
	n_year  <- n_distinct(df$year)
	if (nrow(df) < 6 || n_year < 3 || var(df$conc) == 0) return(NULL)
	df <- mutate(df, year_num = as.numeric(as.character(year)))
	k  <- min(4, n_year)                             # avoid rank-deficiency
	form <- if (n_distinct(df$season) > 1)
		conc ~ season + s(year_num, k = k)
	else  conc ~            s(year_num, k = k)
	tryCatch(gam(form, data = df, method = "REML"), error = function(e) NULL)
}

predict_gam <- function(mod, site, catchment) {
	if (is.null(mod)) return(NULL)
	expand_grid(
		Site      = site,
		Catchment = catchment,
		season    = levels(mod$model$season),
		year_num  = as.numeric(as.character(campaign_years))
	) %>%
		mutate(fit = predict(mod, newdata = tibble(
			season   = season,
			year_num = year_num)),
			year = factor(year_num)) %>%
		select(-year_num)
}

pred_B <- metal_long %>%
	filter(Site %in% hi_res_sites) %>%
	group_by(Site, Catchment, metal) %>%
	nest() %>%
	mutate(model = map(data, safe_gam),
			 pred  = pmap(list(model, Site, Catchment), predict_gam)) %>%
	filter(map_lgl(pred, ~ !is.null(.x))) %>%     # keep only successful fits
	ungroup() %>%                                 # ← key line!
	select(-Site, -Catchment) %>%                 # now they’re really gone
	unnest(pred) %>%                              # no duplicate-name error
	pivot_wider(names_from  = metal,
					values_from = fit,
					names_glue  = "{metal}_pred")

## ── 7  Combine A and B  (B overrides A)───────────────────────────────────────
meta_enriched <- meta %>%
	# 1st: year-specific GAM predictions (hi-res only)
	left_join(pred_B,
				 by = c("Site", "Catchment", "season", "year")) %>%
	# 2nd: fill remaining NA with season means
	left_join(pred_A,
				 by = c("Site", "Catchment", "season"),
				 suffix = c("", ".A")) %>%
	mutate(across(ends_with("_pred"),
					  ~ coalesce(.x, get(paste0(cur_column(), ".A"))))) %>%
	select(-ends_with(".A"))

# Eliminate sites with incomplete metal records
meta_enriched <- meta_enriched %>%
	filter(if_all(ends_with("_pred"), ~ !is.na(.x)))

## ── 8  Example community analysis (PERMANOVA)─────────────────────────────────
vars  <- c("Lead_pred","Zinc_pred","Copper_pred")     # metals you trust
keep  <- complete.cases(meta_enriched[ vars ])
adonis2(
	bcurt ~ Lead_pred + Zinc_pred + Copper_pred + season,
	data   = meta_enriched,
	strata = meta_enriched$Site,
	method = "bray",
	by     = "margin")
