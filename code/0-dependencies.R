if (!require("pacman", quietly = TRUE)) {
	install.packages("pacman")
}

# Packages I want (CRAN first, then Bioconductor)
pkgs <- c(
	"devtools", "vegan", "tidyverse", "report", "lme4", "lmerTest", "permute",
	"DHARMa", "ggforce", "ggpubr", "performance", "AICcPermanova", "furrr",
	"future", "BiocManager", "conflicted"
)
pkgs <- c(
	"devtools", "furrr", "future", "agricolae", "vegan", "report",
	"tidyverse", "permute", "FactoMineR", "factoextra", "lme4",
	"writexl", "viridisLite", "mgcv", "car", "mvabund", "gratia",
	"rstatix", "lubridate", "ggforce", "ggpubr", "performance", "lmerTest",
	"DESeq2", "car", "conflicted", "AICcPermanova", "raster", "DHARMa"
)

# Install (CRAN or Bioconductor) and load in one go
pacman::p_load(char = pkgs, try.bioconductor = TRUE)

BiocManager::version()
BiocManager::install("ANCOMBC")
BiocManager::install("phyloseq")
BiocManager::install("ALDEx2")
BiocManager::install("BiocParallel")
library(BiocParallel)
library(ALDEx2)
library(ANCOMBC)
library(phyloseq)

# set default packages for conflicting function calls. Calls intended for non-default packages must be specified with "::{package_name}" prefix
# conflict_prefer("first", "dplyr")
 conflict_prefer("filter", "dplyr")
# conflict_prefer("rename", "dplyr")
# conflict_prefer("slice", "dplyr")
 conflict_prefer("select", "dplyr")
# conflict_prefer("unique", "base")
# conflict_prefer("as.data.frame", "base")
# conflict_prefer("as.factor", "base")
# conflict_prefer("expand", "tidyr")
# conflict_prefer("getData", "raster")
# conflict_prefer("s", "mgcv")
# conflict_prefer("gam", "mgcv")
# conflicts_prefer(lmerTest::lmer)


# remove if needed
 remove.packages("lmerTest")
