# R script to visualize and statistically analyze data from CoMP-MRS
#
# Authors:
#   Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
#   Diana G. Rotaru, Ph.D. (diana.rotaru@meduniwien.ac.at)
#
# Last updated: 2026-03-30


# Initialize ------------------------------------------------------------------

rm(list = ls()) # Clear the current environment
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE) # If using RStudio, this clears all plots
cat("\014") # CTRL+L (clear console)


# Analysis options ------------------------------------------------------------

show_amcharts <- FALSE # Set to TRUE to create interactive 3D pie charts using amCharts4
show_facet_plots <- FALSE # Set to TRUE to create facet plots of spectral quality metrics by different grouping variables
show_model_diagnostics <- TRUE # Set to TRUE to show model diagnostic plots 
# (e.g., residuals, Q-Q plots) for linear mixed-effects models
save_csv <- TRUE # Set to TRUE to save descriptive statistics tables as CSV files in the derivatives directory


# Set data directories --------------------------------------------------------

##### CHANGE THESE DIRECTORIES AS NEEDED #####

base_dir  <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
data_dir  <- file.path(base_dir, "data")
deriv_dir <- file.path(data_dir, "derivatives")
plots_dir <- file.path(deriv_dir, "plots")

if (!dir.exists(deriv_dir)) {
  dir.create(deriv_dir)
}

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
} else {
  unlink(plots_dir, recursive = TRUE)
  dir.create(plots_dir)
}


# Install/load packages -------------------------------------------------------

# Function to install and/or load necessary packages
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
  "lmerTest"
)

load.packages(packages)


# Set inline functions --------------------------------------------------------

nloptFun <- function(fn, par, lower, upper, control=list(), ...) {
  for (n in names(defaultControl))
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
  res <- nloptr(x0=par, eval_f=fn, lb=lower, ub=upper, opts=control, ...)
  with(res, list(par=solution,
                 fval=objective,
                 feval=iterations,
                 conv=if (status>0) 0 else status,
                 message=message))
}

cs. <- function(x) scale(x, center=T, scale=T) # for standardizing (z-transforming)
# outcome and predictor variables;
# aids with model convergence and
# interpretability of parameter estimates


# Load data -------------------------------------------------------------------

# Load the participants.csv file with appropriate column types
DATA <- read_csv(
  file.path(data_dir, "CoMP_MRS_Rstats_input.csv"),
  col_types = list(
    CompID = col_factor(),
    SiteID = col_factor(),
    DP = col_factor(),
    AnimalID = col_factor(),
    AnimalSpecies = col_factor(),
    AnimalStrain = col_factor(),
    AnimalAge = col_double(),
    AnimalSex = col_factor(),
    AnimalWeight = col_double(),
    MRvendor = col_factor(),
    MRfield = col_double(),
    MRsequence = col_factor(),
    MRbrainregion = col_factor(),
    MRvoxelvolume = col_double(),
    MRaverages = col_double(),
    MRsoftwareversion = col_factor(),
    MRcoildetail = col_factor(),
    MRSsw = col_double(),
    MRSnpts = col_double(),
    MRSTE = col_double(),
    MRSTR = col_double(),
    MRSshimmethod = col_factor(),
    LW = col_double(),
    SNR = col_double(),
    SNR_LW_Ratio = col_double(),
    CompCheck = col_factor()
  )
)

# Remove "compMR" prefix from CompID to make it more concise (optional)
DATA <- DATA %>%
  mutate(CompID = str_remove_all(CompID, "compMR"))


# Spectral quality metrics and normalization ----------------------------------

DATA <- DATA %>%
  mutate(
    SNR_LW_Ratio_norm   = SNR_LW_Ratio / (sqrt(MRaverages) * MRvoxelvolume),
    SNR_LW_Product_norm = (LW * SNR) / (sqrt(MRaverages) * MRvoxelvolume),
    LW_norm             = LW / (sqrt(MRaverages) * MRvoxelvolume),
    SNR_norm            = SNR / (sqrt(MRaverages) * MRvoxelvolume)
  )


# Clean up data ----------------------------------------------------------

# Ensure that all missing values are represented as NA (in case there are any
# other representations of missing data in the original CSV)
DATA[is.na(DATA)] <- NA

### Remove datasets that failed QC --------------------------------------------

DATA <- DATA %>%
  filter(CompCheck != "fail") %>%
  select(-CompCheck) # Remove the CompCheck column as it's no longer needed

### Remove outliers based on the MAD method -----------------------------------
# This will also remove any rows where LW or SNR is NA, which is appropriate for outlier removal

source(file.path(base_dir, "fun", "mad_outlier_removal.R"))

DATA_clean <- mad_outlier_removal(
  data      = DATA,
  outcome   = "SNR_LW_Ratio_norm",
  groups    = "DP",  # outer → inner
  strategy  = "group",  # options: "global", "group", "multilevel"
  threshold = 2.5
)
DATA_orig <- DATA
DATA <- DATA_clean$data_clean # Update DATA to the cleaned dataset

# Optional but recommended:
# aggregate to one row per DP for cleaner plots
DATA_DP <- DATA %>%
  dplyr::group_by(
    DP, SiteID, AnimalSpecies, AnimalSex,
    MRvendor, MRfield, MRsequence, MRbrainregion
  ) %>%
  dplyr::summarise(
    LW_norm = mean(LW_norm, na.rm = TRUE),
    SNR_norm = mean(SNR_norm, na.rm = TRUE),
    SNR_LW_Product_norm = mean(SNR_LW_Product_norm, na.rm = TRUE),
    SNR_LW_Ratio_norm = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
    .groups = "drop"
  )


# Descriptive stats -----------------------------------------------------------

STATS <- list()

### Descriptive stats function --------------------------------------------------

make_stats <- function(data, group_var, value_var, prefix) {
  data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      mean = mean(.data[[value_var]], na.rm = TRUE),
      sd   = sd(.data[[value_var]], na.rm = TRUE),
      cv   = cv(.data[[value_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(
      !!paste0("mean", prefix) := mean,
      !!paste0("sd", prefix)   := sd,
      !!paste0("cv", prefix)   := cv
    )
}

### DP ------------------------------------------------------------------------

STATS$DP <- list(
  LW      = make_stats(DATA, "DP", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "DP", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "DP", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "DP", "SNR_LW_Product_norm", "Product")
)

### SiteID --------------------------------------------------------------------

STATS$SiteID <- list(
  LW      = make_stats(DATA, "SiteID", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "SiteID", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "SiteID", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "SiteID", "SNR_LW_Product_norm", "Product")
)

### AnimalSpecies -------------------------------------------------------------

STATS$AnimalSpecies <- list(
  LW      = make_stats(DATA, "AnimalSpecies", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "AnimalSpecies", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "AnimalSpecies", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "AnimalSpecies", "SNR_LW_Product_norm", "Product")
)

### MRvendor ------------------------------------------------------------------

STATS$MRvendor <- list(
  LW      = make_stats(DATA, "MRvendor", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "MRvendor", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "MRvendor", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "MRvendor", "SNR_LW_Product_norm", "Product")
)

### MRfield -------------------------------------------------------------------

STATS$MRfield <- list(
  LW      = make_stats(DATA, "MRfield", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "MRfield", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "MRfield", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "MRfield", "SNR_LW_Product_norm", "Product")
)

### MRsequence ----------------------------------------------------------------

STATS$MRsequence <- list(
  LW      = make_stats(DATA, "MRsequence", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "MRsequence", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "MRsequence", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "MRsequence", "SNR_LW_Product_norm", "Product")
)

### MRbrainregion -------------------------------------------------------------

STATS$MRbrainregion <- list(
  LW      = make_stats(DATA, "MRbrainregion", "LW_norm", "LW"),
  SNR     = make_stats(DATA, "MRbrainregion", "SNR_norm", "SNR"),
  Ratio   = make_stats(DATA, "MRbrainregion", "SNR_LW_Ratio_norm", "Ratio"),
  Product = make_stats(DATA, "MRbrainregion", "SNR_LW_Product_norm", "Product")
)

### AnimalSpecies x MRvendor --------------------------------------------------

STATS$AnimalSpecies_MRvendor <- DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(
    meanLW      = mean(LW_norm, na.rm = TRUE),
    sdLW        = sd(LW_norm, na.rm = TRUE),
    cvLW        = cv(LW_norm, na.rm = TRUE),
    meanSNR     = mean(SNR_norm, na.rm = TRUE),
    sdSNR       = sd(SNR_norm, na.rm = TRUE),
    cvSNR       = cv(SNR_norm, na.rm = TRUE),
    meanRatio   = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
    sdRatio     = sd(SNR_LW_Ratio_norm, na.rm = TRUE),
    cvRatio     = cv(SNR_LW_Ratio_norm, na.rm = TRUE),
    meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE),
    sdProduct   = sd(SNR_LW_Product_norm, na.rm = TRUE),
    cvProduct   = cv(SNR_LW_Product_norm, na.rm = TRUE),
    .groups = "drop"
  )

### AnimalSpecies x MRbrainregion ---------------------------------------------

STATS$AnimalSpecies_MRbrainregion <- DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(
    meanLW      = mean(LW_norm, na.rm = TRUE),
    sdLW        = sd(LW_norm, na.rm = TRUE),
    cvLW        = cv(LW_norm, na.rm = TRUE),
    meanSNR     = mean(SNR_norm, na.rm = TRUE),
    sdSNR       = sd(SNR_norm, na.rm = TRUE),
    cvSNR       = cv(SNR_norm, na.rm = TRUE),
    meanRatio   = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
    sdRatio     = sd(SNR_LW_Ratio_norm, na.rm = TRUE),
    cvRatio     = cv(SNR_LW_Ratio_norm, na.rm = TRUE),
    meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE),
    sdProduct   = sd(SNR_LW_Product_norm, na.rm = TRUE),
    cvProduct   = cv(SNR_LW_Product_norm, na.rm = TRUE),
    .groups = "drop"
  )

### All -----------------------------------------------------------------------

STATS$all <- list(
  meanLW      = mean(DATA$LW_norm, na.rm = TRUE),
  sdLW        = sd(DATA$LW_norm, na.rm = TRUE),
  cvLW        = cv(DATA$LW_norm, na.rm = TRUE),
  meanSNR     = mean(DATA$SNR_norm, na.rm = TRUE),
  sdSNR       = sd(DATA$SNR_norm, na.rm = TRUE),
  cvSNR       = cv(DATA$SNR_norm, na.rm = TRUE),
  meanRatio   = mean(DATA$SNR_LW_Ratio_norm, na.rm = TRUE),
  sdRatio     = sd(DATA$SNR_LW_Ratio_norm, na.rm = TRUE),
  cvRatio     = cv(DATA$SNR_LW_Ratio_norm, na.rm = TRUE),
  meanProduct = mean(DATA$SNR_LW_Product_norm, na.rm = TRUE),
  sdProduct   = sd(DATA$SNR_LW_Product_norm, na.rm = TRUE),
  cvProduct   = cv(DATA$SNR_LW_Product_norm, na.rm = TRUE)
)


# Save descriptive stats tables -----------------------------------------------

if (save_csv) {
  
  write_csv(STATS$DP$LW, file.path(deriv_dir, "stats_DP_LW.csv"))
  write_csv(STATS$DP$SNR, file.path(deriv_dir, "stats_DP_SNR.csv"))
  write_csv(STATS$DP$Ratio, file.path(deriv_dir, "stats_DP_Ratio.csv"))
  write_csv(STATS$DP$Product, file.path(deriv_dir, "stats_DP_Product.csv"))
  
  write_csv(STATS$SiteID$LW, file.path(deriv_dir, "stats_SiteID_LW.csv"))
  write_csv(STATS$SiteID$SNR, file.path(deriv_dir, "stats_SiteID_SNR.csv"))
  write_csv(STATS$SiteID$Ratio, file.path(deriv_dir, "stats_SiteID_Ratio.csv"))
  write_csv(STATS$SiteID$Product, file.path(deriv_dir, "stats_SiteID_Product.csv"))
  
  write_csv(STATS$AnimalSpecies$LW, file.path(deriv_dir, "stats_AnimalSpecies_LW.csv"))
  write_csv(STATS$AnimalSpecies$SNR, file.path(deriv_dir, "stats_AnimalSpecies_SNR.csv"))
  write_csv(STATS$AnimalSpecies$Ratio, file.path(deriv_dir, "stats_AnimalSpecies_Ratio.csv"))
  write_csv(STATS$AnimalSpecies$Product, file.path(deriv_dir, "stats_AnimalSpecies_Product.csv"))
  
  write_csv(STATS$MRvendor$LW, file.path(deriv_dir, "stats_MRvendor_LW.csv"))
  write_csv(STATS$MRvendor$SNR, file.path(deriv_dir, "stats_MRvendor_SNR.csv"))
  write_csv(STATS$MRvendor$Ratio, file.path(deriv_dir, "stats_MRvendor_Ratio.csv"))
  write_csv(STATS$MRvendor$Product, file.path(deriv_dir, "stats_MRvendor_Product.csv"))
  
  write_csv(STATS$MRfield$LW, file.path(deriv_dir, "stats_MRfield_LW.csv"))
  write_csv(STATS$MRfield$SNR, file.path(deriv_dir, "stats_MRfield_SNR.csv"))
  write_csv(STATS$MRfield$Ratio, file.path(deriv_dir, "stats_MRfield_Ratio.csv"))
  write_csv(STATS$MRfield$Product, file.path(deriv_dir, "stats_MRfield_Product.csv"))
  
  write_csv(STATS$MRsequence$LW, file.path(deriv_dir, "stats_MRsequence_LW.csv"))
  write_csv(STATS$MRsequence$SNR, file.path(deriv_dir, "stats_MRsequence_SNR.csv"))
  write_csv(STATS$MRsequence$Ratio, file.path(deriv_dir, "stats_MRsequence_Ratio.csv"))
  write_csv(STATS$MRsequence$Product, file.path(deriv_dir, "stats_MRsequence_Product.csv"))
  
  write_csv(STATS$MRbrainregion$LW, file.path(deriv_dir, "stats_MRbrainregion_LW.csv"))
  write_csv(STATS$MRbrainregion$SNR, file.path(deriv_dir, "stats_MRbrainregion_SNR.csv"))
  write_csv(STATS$MRbrainregion$Ratio, file.path(deriv_dir, "stats_MRbrainregion_Ratio.csv"))
  write_csv(STATS$MRbrainregion$Product, file.path(deriv_dir, "stats_MRbrainregion_Product.csv"))
  
  write_csv(
    STATS$AnimalSpecies_MRvendor,
    file.path(deriv_dir, "stats_AnimalSpecies_MRvendor.csv")
  )
  
  write_csv(
    STATS$AnimalSpecies_MRbrainregion,
    file.path(deriv_dir, "stats_AnimalSpecies_MRbrainregion.csv")
  )
  
  write_csv(
    tibble(
      meanLW = STATS$all$meanLW,
      sdLW = STATS$all$sdLW,
      cvLW = STATS$all$cvLW,
      meanSNR = STATS$all$meanSNR,
      sdSNR = STATS$all$sdSNR,
      cvSNR = STATS$all$cvSNR,
      meanRatio = STATS$all$meanRatio,
      sdRatio = STATS$all$sdRatio,
      cvRatio = STATS$all$cvRatio,
      meanProduct = STATS$all$meanProduct,
      sdProduct = STATS$all$sdProduct,
      cvProduct = STATS$all$cvProduct
    ),
    file.path(deriv_dir, "stats_all.csv")
  )
  
}


# Pie chart function ----------------------------------------------------------

pie_dir <- file.path(plots_dir, "piecharts")
dir.create(pie_dir, recursive = TRUE, showWarnings = FALSE)

# Custom color palette 
pie_palette <- c(
  "#4C77C2",  # blue
  "#EA7E2D",  # orange
  "#A7A7A7",  # gray
  "#F2BE00",  # yellow
  "#5B9BD5",  # light blue
  "#70AD47",  # green
  "#264478",  # dark blue
  "#9E480E",  # brown
  "#636363",  # dark gray
  "#997300",  # olive
  "#255E91",  # steel blue
  "#43682B"   # darker green
)

make_pie_chart <- function(data, var, plot_title = NULL, file_name = NULL, output_dir) {
  var_sym <- rlang::sym(var)
  
  # AnimalSex should be based on row count, not unique DP-SiteID combinations
  if (var == "AnimalSex") {
    plot_data <- data %>%
      dplyr::filter(!is.na(!!var_sym)) %>%
      dplyr::group_by(!!var_sym) %>%
      dplyr::summarise(n_count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        perc = n_count / sum(n_count),
        perc_label = scales::percent(perc, accuracy = 0.1),
        legend_label = paste0(as.character(!!var_sym), " (n=", n_count, ", ", perc_label, ")")
      ) %>%
      dplyr::arrange(dplyr::desc(n_count))
  } else {
    # All other pie charts should be based on unique DP-SiteID combinations
    plot_data <- data %>%
      dplyr::filter(!is.na(!!var_sym), !is.na(DP), !is.na(SiteID)) %>%
      dplyr::distinct(DP, SiteID, !!var_sym) %>%
      dplyr::group_by(!!var_sym) %>%
      dplyr::summarise(n_count = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(
        perc = n_count / sum(n_count),
        perc_label = scales::percent(perc, accuracy = 0.1),
        legend_label = paste0(as.character(!!var_sym), " (n=", n_count, ", ", perc_label, ")")
      ) %>%
      dplyr::arrange(dplyr::desc(n_count))
  }
  
  if (is.null(plot_title)) plot_title <- paste(var, "distribution")
  if (is.null(file_name)) file_name <- paste0(var, "_piechart.png")
  
  fill_vals <- rep(pie_palette, length.out = nrow(plot_data))
  
  p <- ggplot(plot_data, aes(x = 1, y = perc, fill = legend_label)) +
    geom_col(
      color = "white",
      linewidth = 0.4,
      width = 1
    ) +
    geom_text(
      aes(x = 1.3, label = perc_label), # aes(label = count_label),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 5,
      fontface = "bold",
      family = "serif"
    ) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = fill_vals) +
    xlim(0.5, 1.5) +
    theme_void() +
    labs(
      title = plot_title,
      fill = NULL
    ) +
    theme(
      plot.title = element_text(
        family = "serif",
        face = "plain",
        size = 18,
        hjust = 0.5,
        color = "#4D4D4D"
      ),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(
        family = "serif",
        size = 10,
        color = "#4D4D4D"
      ),
      legend.key.size = unit(0.35, "cm"),
      legend.spacing.x = unit(0.2, "cm"),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    guides(
      fill = guide_legend(
        nrow = ceiling(nrow(plot_data) / 2),
        byrow = TRUE
      )
    )
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, file_name),
    plot = p,
    height = 5.8,
    width = 5.8,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  return(p)
}


# Create pie charts with percentages ------------------------------------------

pie_specs <- list(
  list(var = "SiteID",         title = "Site (#DPs)"),
  list(var = "AnimalSpecies",  title = "Species (#DPs)"),
  list(var = "AnimalStrain",   title = "Animal strain (#DPs)"),
  list(var = "AnimalSex",      title = "Sex (#animals)"),
  list(var = "MRvendor",       title = "Vendor (#DPs)"),
  list(var = "MRfield",        title = "Field strength (#DPs)"),
  list(var = "MRsequence",     title = "Pulse sequence (#DPs)"),
  list(var = "MRbrainregion",  title = "Voxel location (#DPs)"),
  list(var = "MRaverages",     title = "Nr. of averages (#DPs)")
)

pie_plots <- lapply(pie_specs, function(x) {
  make_pie_chart(
    data = DATA,
    var = x$var,
    plot_title = x$title,
    file_name = paste0(x$var, "_piechart.png"),
    output_dir = pie_dir
  )
})

names(pie_plots) <- sapply(pie_specs, `[[`, "var")


# Combined preview of pie charts ----------------------------------------------

combined_pies <- patchwork::wrap_plots(pie_plots, ncol = 3)
ggsave(
  filename = file.path(plots_dir, "all_piecharts_combined.png"),
  plot = combined_pies,
  width = 16,
  height = 18,
  units = "in",
  dpi = 300
)


# 3D pie chart function -------------------------------------------------------

# Source - https://stackoverflow.com/a/76142219
# Posted by Stéphane Laurent
# Retrieved 2026-03-26, License - CC BY-SA 4.0

make_pie_data_amcharts <- function(data, var) {
  var_sym <- rlang::sym(var)
  
  if (var == "AnimalSex") {
    plot_data <- data %>%
      dplyr::filter(!is.na(!!var_sym)) %>%
      dplyr::group_by(!!var_sym) %>%
      dplyr::summarise(value = dplyr::n(), .groups = "drop")
  } else {
    plot_data <- data %>%
      dplyr::filter(!is.na(!!var_sym), !is.na(DP), !is.na(SiteID)) %>%
      dplyr::distinct(DP, SiteID, !!var_sym) %>%
      dplyr::group_by(!!var_sym) %>%
      dplyr::summarise(value = dplyr::n(), .groups = "drop")
  }
  
  plot_data %>%
    dplyr::rename(group = !!var_sym) %>%
    dplyr::mutate(group = as.character(group))
}


# Create 3D pie charts with percentages --------------------------------------

if (show_amcharts) {
  
  amcharts_dir <- file.path(plots_dir, "amcharts_pies")
  dir.create(amcharts_dir, recursive = TRUE, showWarnings = FALSE)
  
  vars <- c(
    "SiteID",
    "AnimalSpecies",
    "AnimalStrain",
    "AnimalSex",
    "MRvendor",
    "MRfield",
    "MRsequence",
    "MRbrainregion",
    "MRaverages"
  )
  
  charts <- list()
  
  for (v in vars) {
    dat <- make_pie_data_amcharts(DATA, v) %>%
      dplyr::arrange(dplyr::desc(value))
    
    charts[[v]] <- amPieChart(
      data = dat,
      category = "group",
      value = "value",
      threeD = TRUE,
      variableDepth = TRUE,
      chartTitle = v,
      legend = TRUE
    )
    
    charts[[v]] <- htmlwidgets::onRender(
      charts[[v]],
      "
    function(el, x) {
      var chart = this.chart;

      if (chart.series && chart.series.values.length > 0) {
        var series = chart.series.values[0];

        series.startAngle = 300;
        series.endAngle   = 660;

        series.slices.template.stroke = am4core.color('#FFFFFF');
        series.slices.template.strokeWidth = 2;
        series.slices.template.strokeOpacity = 1;
      }
    }
    "
    )
    
    print(charts[[v]])
    
    htmlwidgets::saveWidget(
      widget = charts[[v]],
      file = file.path(amcharts_dir, paste0(v, "_piechart.html")),
      selfcontained = TRUE
    )
    
    browseURL(file.path(amcharts_dir,  paste0(v, "_piechart.html")))
  }
  
}


# Plot options -----------------------------------------------------------

# Common theme
theme_comp <- function() {
  theme_classic2() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        color = "black",
        family = "Arial"
      ),
      axis.title = element_text(
        face = "bold",
        size = 13,
        color = "black",
        family = "Arial"
      ),
      axis.text = element_text(
        size = 11,
        color = "black",
        family = "Arial"
      ),
      strip.text = element_text(
        face = "bold",
        size = 12,
        color = "black",
        family = "Arial"
      ),
      legend.position = "none",
      axis.line = element_line(
        color = "black",
        linewidth = 0.4,
        lineend = "square"
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.4
      ),
      panel.spacing = unit(1, "lines")
    )
}

# Dynamic color function for any number of groups
get_discrete_colors <- function(n) {
  if (n == 1) {
    return(c("#4C77C2")) # blue
  } else if (n == 2) {
    return(c("#4C77C2", "#EA7E2D")) # blue and orange
  } else if (n > 2 && n <= 8) {
    return(RColorBrewer::brewer.pal(n, "Set2"))
  } else {
    return(grDevices::hcl.colors(n, palette = "Dynamic"))
  }
}


# Dot plot function ------------------------------------------------------

make_dot_plot <- function(data, x_var, y_var, y_label, group_var, file_name) {
  
  every_nth <- function(n) {
    function(x) x[seq(1, length(x), by = n)]
  }
  
  data <- data %>%
    arrange(.data[[group_var]], .data[[x_var]]) %>%
    mutate(!!x_var := factor(.data[[x_var]], levels = unique(.data[[x_var]])))
  
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], color = factor(.data[[group_var]]))) +
    geom_point(
      size = 1
    ) +
    labs(
      title = paste(y_label, "by", x_var),
      x = x_var,
      y = y_label,
      color = group_var
    ) +
    scale_x_discrete(
      breaks = every_nth(10) # Show every 10th label on x-axis
    ) +
    theme_comp() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  ggsave(
    filename = file.path(plots_dir, file_name),
    plot = p,
    width = 7,
    height = 5,
    units = "in",
    dpi = 300
  )
  
  return(p)
}


# Violin plots function -------------------------------------------------------

# General violin plot function
make_violin_plot <- function(data, x_var, y_var, y_label, file_name, level_subset = NULL) {
  
  plot_data <- data %>%
    dplyr::filter(!is.na(.data[[x_var]]), !is.na(.data[[y_var]])) %>%
    dplyr::mutate(
      x_group = as.factor(.data[[x_var]])
    )
  
  # Optional subset of levels for crowded variables
  if (!is.null(level_subset)) {
    plot_data <- plot_data %>%
      dplyr::filter(x_group %in% level_subset) %>%
      dplyr::mutate(x_group = factor(x_group, levels = level_subset))
  }
  
  n_groups <- nlevels(plot_data$x_group)
  fill_colors <- get_discrete_colors(n_groups)
  
  p <- ggplot(plot_data, aes(x = x_group, y = .data[[y_var]], fill = x_group)) +
    # geom_violin(
    #   position = position_dodge(0),
    #   width = 1,
    #   trim = FALSE,
    #   alpha = 0.8,
    #   linewidth = 0.25,
    #   color = "black") +
    geom_boxplot(
      width = 0.65,
      # fill = "#EA7E2D",
      alpha = 0.9,
      outlier.shape = 21,
      outlier.size = 1.8,
      outlier.stroke = 0.3,
      color = "black",
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.3,
      size = 1.2,
      color = "black"
    ) +
    geom_point(
      size = 1
    ) +
    labs(
      title = paste(y_label, "by", x_var),
      x = x_var,
      y = y_label
    ) +
    scale_fill_manual(
      values = fill_colors,
      na.translate = FALSE
    ) +
    theme_comp()
  
  # Rotate labels for crowded x-axes
  if (nlevels(plot_data$x_group) > 4) {
    p <- p + theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  }
  
  # Make wider plots for many categories
  plot_width <- if (nlevels(plot_data$x_group) > 10) 10 else 7
  
  print(p)
  
  ggsave(
    filename = file.path(plots_dir, file_name),
    plot = p,
    width = plot_width,
    height = 5,
    units = "in",
    dpi = 300
  )
  
  return(p)
}


# Create plots ----------------------------------------------------------------

### Dot plots -----------------------------------------------------------------

dotplot_orig <- make_dot_plot(
  data = DATA_orig,   # use DATA here instead if you want row-level plots
  x_var = "CompID",
  y_var = "SNR_LW_Ratio_norm",
  y_label = "Normalized SNR/LW ratio",
  group_var = "DP",
  file_name = "dotplot_SNR_LW_Ratio_norm_by_Subj_orig.png"
)

dotplot <- make_dot_plot(
  data = DATA,   # use DATA here instead if you want row-level plots
  x_var = "CompID",
  y_var = "SNR_LW_Ratio_norm",
  y_label = "Normalized SNR/LW ratio",
  group_var = "DP",
  file_name = "dotplot_SNR_LW_Ratio_norm_by_Subj.png"
)

### Y variables (measured) ----------------------------------------------------

y_specs <- list(
  # list(var = "LW_norm",             label = "Normalized LW"),
  # list(var = "SNR_norm",            label = "Normalized SNR"),
  # list(var = "SNR_LW_Product_norm", label = "Normalized SNR×LW"),
  list(var = "SNR_LW_Ratio_norm",   label = "Normalized SNR/LW")
)

### X variables (grouping factors) --------------------------------------------

x_vars <- c(
  "DP",
  "SiteID",
  "AnimalSpecies",
  "AnimalSex",
  "MRvendor",
  "MRfield",
  "MRsequence",
  "MRbrainregion"
)

### Violin plots --------------------------------------------------------------

all_violin_plots <- list()

for (y in y_specs) {
  for (x in x_vars) {
    
    # Special handling for crowded variables: split into two plots
    if (x %in% c("DP")) {
      
      all_levels <- DATA_DP %>%
        dplyr::filter(!is.na(.data[[x]])) %>%
        dplyr::pull(.data[[x]]) %>%
        unique() %>%
        as.character() %>%
        sort()
      
      split_index <- ceiling(length(all_levels) / 2)
      levels_part1 <- all_levels[1:split_index]
      levels_part2 <- all_levels[(split_index + 1):length(all_levels)]
      
      # First half
      file_name_1 <- paste0("violinplot_", y$var, "_by_", x, "_part1.png")
      key_1 <- paste(y$var, x, "part1", sep = "_by_")
      
      all_violin_plots[[key_1]] <- make_violin_plot(
        data = DATA_DP,
        x_var = x,
        y_var = y$var,
        y_label = paste0(y$label, " (", x, " 1)"),
        file_name = file_name_1,
        level_subset = levels_part1
      )
      
      # Second half
      if (length(levels_part2) > 0) {
        file_name_2 <- paste0("violinplot_", y$var, "_by_", x, "_part2.png")
        key_2 <- paste(y$var, x, "part2", sep = "_by_")
        
        all_violin_plots[[key_2]] <- make_violin_plot(
          data = DATA_DP,
          x_var = x,
          y_var = y$var,
          y_label = paste0(y$label, " (", x, " 2)"),
          file_name = file_name_2,
          level_subset = levels_part2
        )
      }
      
    } else {
      
      file_name <- paste0("violinplot_", y$var, "_by_", x, ".png")
      key <- paste(y$var, x, sep = "_by_")
      
      all_violin_plots[[key]] <- make_violin_plot(
        data = DATA_DP,
        x_var = x,
        y_var = y$var,
        y_label = y$label,
        file_name = file_name
      )
    }
  }
}

### Double variable boxplots: MRvendor and AnimalSpecies ----------------------

# This will create boxplots of each y variable, faceted by MRvendor and AnimalSpecies
# You can modify the function to facet by other combinations of variables if desired

make_boxplot_vendor_species <- function(data, y_var, y_label, file_name) {
  
  plot_data <- data %>%
    dplyr::filter(
      !is.na(MRvendor),
      !is.na(AnimalSpecies),
      !is.na(.data[[y_var]])
    ) %>%
    dplyr::mutate(
      MRvendor = as.factor(MRvendor),
      AnimalSpecies = as.factor(AnimalSpecies)
    )
  
  p <- ggplot(plot_data, aes(x = MRvendor, y = .data[[y_var]], fill = MRvendor)) +
    geom_violin(
      position = position_dodge(0),
      width = 1,
      trim = FALSE,
      alpha = 0.8,
      linewidth = 0.25,
      color = "black") +
    geom_boxplot(
      width = 0.65,
      alpha = 0.9,
      outlier.shape = 21,
      outlier.size = 2.0,
      outlier.stroke = 0.3,
      color = "black",
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.12,
      alpha = 0.4,
      size = 1.6,
      color = "black"
    ) +
    geom_point(
      size = 1
    ) +
    facet_wrap(~ AnimalSpecies, scales = "free_x") +
    labs(
      title = paste(y_label, "by vendor and animal species"),
      x = "Vendor",
      y = y_label
    ) +
    scale_fill_brewer(palette = "Set2", na.translate = FALSE) +
    theme_comp()
  
  print(p)
  
  ggsave(
    filename = file.path(plots_dir, file_name),
    plot = p,
    width = 8,
    height = 5,
    units = "in",
    dpi = 300
  )
  
  return(p)
}

if (show_facet_plots) {
  
  facet_specs <- list(
    list(
      var = "LW_norm",
      label = "Normalized LW",
      file = "boxplot_LW_norm_by_MRvendor_AnimalSpecies.png"
    ),
    list(
      var = "SNR_norm",
      label = "Normalized SNR",
      file = "boxplot_SNR_norm_by_MRvendor_AnimalSpecies.png"
    ),
    list(
      var = "SNR_LW_Product_norm",
      label = "Normalized SNR×LW product",
      file = "boxplot_SNR_LW_Product_norm_by_MRvendor_AnimalSpecies.png"
    ),
    list(
      var = "SNR_LW_Ratio_norm",
      label = "Normalized SNR/LW ratio",
      file = "boxplot_SNR_LW_Ratio_norm_by_MRvendor_AnimalSpecies.png"
    )
  )
  
  facet_plots <- lapply(facet_specs, function(x) {
    make_boxplot_vendor_species(
      data = DATA_DP,   # use DATA here instead if you want row-level plots
      y_var = x$var,
      y_label = x$label,
      file_name = x$file
    )
  })
  
  names(facet_plots) <- sapply(facet_specs, `[[`, "var")
  
}


# SNR_LW_ratio_norm -----------------------------------------------------------

### Linear-mixed models -------------------------------------------------------

# Selection of optimizers to use
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA", xtol_rel=1e-6, maxeval=1e5)
optWrap.nloptwrap.NLOPT_LN_BOBYQA <- lmerControl(optimizer="nloptFun")
optWrap.optimx.nlminb             <- lmerControl(optimizer="optimx", optCtrl=list(method="nlminb", eval.max=1e5))
optWrap.minqa.bobyqa              <- lmerControl(optimizer="bobyqa")

optimToUse <- optWrap.nloptwrap.NLOPT_LN_BOBYQA

# Remove rows with missing values in the variables of interest for the models
DATA <- drop_na(DATA, "SNR_LW_Ratio_norm", "MRfield", "DP")

M.SNRLWrationorm.0.0 <- lm(SNR_LW_Ratio_norm ~ 1, data = DATA)
M.SNRLWrationorm.0.1 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | DP),
                             data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.0.1.a <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | MRbrainregion),
                                   data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.0.1.b <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | MRsequence),
                                     data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.0.2 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor), data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.0.3 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | SiteID), data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.0.4 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | DP) + (1 | SiteID), data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.0.full <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | DP),     # (1 | MRvendor) + (1 | SiteID) only account for negligible variance
                                  data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.1.0 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | AnimalSpecies) + (1 | DP), # (1 | AnimalSpecies) only accounts for negligible variance
                             data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.2.0 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | MRbrainregion) + (1 | DP), # (1 | MRbrainregion) only accounts for negligible variance
                                  data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.3.0 <- lme4::lmer(SNR_LW_Ratio_norm ~ (1 | MRsequence) + (1 | DP),
                                  data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.4.0 <- lme4::lmer(SNR_LW_Ratio_norm ~ MRfield + (1 | DP),
                               data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.4.1 <- lme4::lmer(SNR_LW_Ratio_norm ~ MRfield + (MRfield | DP),
                                  data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.5.0 <- lmer(SNR_LW_Ratio_norm ~  AnimalAge + (1 | DP),
                                  data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.5.1 <- lmer(SNR_LW_Ratio_norm ~  AnimalAge + (AnimalAge | DP),
                             data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.6.0 <- lmer(SNR_LW_Ratio_norm ~  AnimalSex + (1 | DP),
                                  data = DATA, REML = FALSE, control = optimToUse)
M.SNRLWrationorm.6.1 <- lmer(SNR_LW_Ratio_norm ~  AnimalSex + (AnimalSex | DP),
                             data = DATA, REML = FALSE, control = optimToUse)

### Model diagnostics ---------------------------------------------------------

model <- M.SNRLWrationorm.3.0

if (show_model_diagnostics) {
  
  cat("\n── Random-effects diagnostics ──\n")
  print(dotplot(ranef(model, condVar = TRUE), scales = "free"))
  
  cat("\n── Model performance indices ──\n")
  print(model_performance(model))
  
  cat("\n── Intraclass Correlation Coefficient (ICC) ──\n")
  print(icc(model))
  
  cat("\n── Check model assumptions (performance) ──\n")
  chk <- check_model(model)
  print(chk)
  
}

### Variance components -------------------------------------------------------

# Save variance components of three-level null model
# print(VarCorr(M.SNRLWrationorm.0.1), comp = c("Variance", "Std.Dev."))
vc <- as.data.frame(VarCorr(M.SNRLWrationorm.0.full))[4]
vpc.M.0.full <- round(vc/sum(vc)*100,1)
dimnames(vpc.M.0.full)[[1]] <- c("DP", "Residual")

# print(VarCorr(M.SNRLWrationorm.0.1), comp = c("Variance", "Std.Dev."))
vc <- as.data.frame(VarCorr(M.SNRLWrationorm.3.0))[4]
vpc.M.3.0 <- round(vc/sum(vc)*100,1)
dimnames(vpc.M.3.0)[[1]] <- c("MRsequence","DP", "Residual")


### Inferential tests ---------------------------------------------------------
# Compare null model to model with DP random effect
LRT.1 <- anova(M.SNRLWrationorm.0.full, M.SNRLWrationorm.3.0)

# Bootstrapping confidence intervals for fixed effects (intercept)
# b.par1 <- bootMer(Y.M0.3, fixef, nsim=1e4) # bootstrap fixed effects
# b.par2 <- bootMer(Y.M1.3, fixef, nsim=1e4)
# boot.ci(b.par1, conf = 0.95, type = "bca", index = 1)
# boot.ci(b.par2, conf = 0.95, type = "bca", index = 1)

confint(M.SNRLWrationorm.4.0, parm = c(3,4), level = 0.95, method = "boot",
        nsim = 500, boot.type = "perc") # alternative bootstrapping approach

# Set up parallel computation for bootstrapping
nc <- detectCores()
clus <- makeCluster(rep("localhost", nc))
clusterEvalQ(clus, {
  library(lme4)
  library(pbkrtest)
})
clusterExport(clus, c("nloptr", "defaultControl", "nloptFun",
                      "DATA", "M.SNRLWrationorm.0.full", "M.SNRLWrationorm.3.0"))

# Bootstrap likelihood ratio tests
PB_LRT.1 <- PBmodcomp(M.SNRLWrationorm.3.0, M.SNRLWrationorm.0.full, nsim = 2e3, cl = clus)

stopCluster(clus)

# End -------------------------------------------------------------------------
