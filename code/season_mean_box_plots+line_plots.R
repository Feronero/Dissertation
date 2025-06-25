## ── 0  Libraries ──────────────────────────────────────────────────────────────
library(tidyverse)      # ggplot2, dplyr, tidyr, purrr
library(lubridate)      # date helpers

## ── 1  Read data ─────────────────────────────────────────────────────────────
metal_raw <- metal_wide    # adjust path if needed

## ── 2  Tidy to long format & add helper columns ──────────────────────────────
metal_long <- metal_raw %>% 
	pivot_longer(Lead:Calcium, names_to = "metal", values_to = "conc",
					 values_drop_na = TRUE) %>% 
	mutate(
		Site      = factor(Site),
		season    = factor(Season, levels = c("Winter","Spring","Summer","Autumn")),
		# build a mid-month Date for continuous‐time plots; change if you have day info
		Date      = make_date(year, month, 15)
	)

# All metals present in the file (16 expected)
metal_levels <- levels(factor(metal_long$metal))

# Ensure facets always appear in the same order
metal_long <- metal_long %>% mutate(metal = factor(metal, levels = metal_levels))

## ── 3  Functions to generate & save plots per site ───────────────────────────
plot_box_season  <- function(one_site_df, out_path) {
	ggplot(one_site_df, aes(x = season, y = conc, group = season)) +
		geom_boxplot(outlier.alpha = .3) +
		facet_wrap(~ metal, scales = "free_y", ncol = 4) +
		labs(title  = paste("Seasonal metal variability – Site", unique(one_site_df$Site)),
			  x = "Season", y = "Concentration") +
		theme_bw() +
		theme(strip.text = element_text(size = 7)) -> p
	ggsave(out_path, p, width = 210, height = 297, units = "mm")  # A4 portrait
}

plot_line_series <- function(one_site_df, out_path) {
	ggplot(one_site_df, aes(Date, conc)) +
		geom_line(na.rm = TRUE) +
		facet_wrap(~ metal, scales = "free_y", ncol = 4) +
		labs(title = paste("Metal time-series – Site", unique(one_site_df$Site)),
			  x = "Date", y = "Concentration") +
		theme_bw() +
		theme(strip.text = element_text(size = 7)) -> p
	ggsave(out_path, p, width = 297, height = 210, units = "mm")  # A4 landscape
}

## ── 4  Loop over sites and write PDFs ────────────────────────────────────────
dir.create("plots_box",  showWarnings = FALSE)
dir.create("plots_line", showWarnings = FALSE)

metal_long %>% 
	group_split(Site) %>%                         # 22 tibbles, one per site
	walk(function(df_site) {
		site_id <- sprintf("%02d", as.integer(df_site$Site[1]))
		
		plot_box_season( df_site,
							  file.path("plots_box",
							  			 paste0("Site_", site_id, "_season_box.pdf")))
		
		plot_line_series(df_site,
							  file.path("plots_line",
							  			 paste0("Site_", site_id, "_time_series.pdf")))
	})