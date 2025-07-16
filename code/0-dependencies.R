if (!require("pacman", quietly = TRUE)) {
	install.packages("pacman")
}

# Packages you want (CRAN first, then Bioconductor)
pkgs <- c(
	"furrr", "future", "agricolae", "vegan",
	"tidyverse", "permute", "FactoMineR", "factoextra",
	"writexl", "viridisLite", "mgcv", "car",
	"rstatix", "lubridate", "ggforce", "ggpubr",
	"DESeq2", "car", "conflicted", "AICcPermanova"
)

# Install (CRAN or Bioconductor) and load in one go
pacman::p_load(char = pkgs, try.bioconductor = TRUE)

# set default packages for conflicting function calls. Calls intended for non-default packages must be specified with "::{package_name}" prefix
conflict_prefer("first", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("slice", "dplyr")
