# Function to install and/or load necessary packages for the analysis.
# It checks if each package is already installed, and if not, it installs it
# before loading it into the R session.

load.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Packages needed for this script
packages <- c(
  "tidyverse",
  "psych",
  "MASS",
  "lme4",
  "emmeans",
  "multcomp",
  "lattice",
  "car",
  "berryFunctions",
  "patchwork",
  "cowplot",
  "ggpubr",
  "irr",
  "DT",
  "raster",
  "rstatix",
  "extrafont",
  "hrbrthemes",
  "see",
  "latex2exp",
  "boot",
  "pbkrtest",
  "parallel",
  "nloptr",
  "optimx",
  "sjPlot",
  "performance",
  "rAmCharts4",
  "htmlwidgets",
  "lmerTest",
  "modelsummary",
  "flextable",
  "officer",
  "svglite"
)

load.packages(packages)
