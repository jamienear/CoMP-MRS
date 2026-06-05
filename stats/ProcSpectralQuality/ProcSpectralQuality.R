# Wrapper R script to analyze data from CoMP-MRS
#
# Authors:
#   Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
#   Diana G. Rotaru, Ph.D. (diana.rotaru@meduniwien.ac.at)
#
# AI Disclosure: Part of this code was generated using Claude Sonnet 4.6;
# the authors then modified and extended the code as needed for the specific
# analyses and visualizations in this project.
#
# Last updated: 2026-04-16

# Initialize ------------------------------------------------------------------

rm(list = ls()) # Clear the current environment
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE) # If using RStudio, this clears all plots
cat("\014") # CTRL+L (clear console)


# Set data directories --------------------------------------------------------

##### CHANGE THESE DIRECTORIES AS NEEDED #####

base_dir  <- dirname(rstudioapi::getSourceEditorContext()$path)
data_dir  <- file.path(base_dir, "data")
deriv_dir <- file.path(data_dir, "derivatives")
plots_dir <- file.path(deriv_dir, "plots")
csv_dir   <- file.path(deriv_dir, "csv")

if (!dir.exists(deriv_dir)) {
  dir.create(deriv_dir)
}

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
} else {
  unlink(plots_dir, recursive = TRUE)
  dir.create(plots_dir)
}

if (!dir.exists(csv_dir)) {
  dir.create(csv_dir)
} else {
  unlink(csv_dir, recursive = TRUE)
  dir.create(csv_dir)
}


# Pre-load code ---------------------------------------------------------------

source(file.path(base_dir, "fun", "mad_outlier_removal.R"))
source(file.path(base_dir, "scripts", "LoadPackages.R"))
source(file.path(base_dir, "scripts", "InlineFun.R"))
source(file.path(base_dir, "scripts", "LoadData.R"))
source(file.path(base_dir, "scripts", "DescripStats.R"))
source(file.path(base_dir, "scripts", "SaveCSV.R"))
source(file.path(base_dir, "scripts", "PlotPieCharts.R"))
source(file.path(base_dir, "scripts", "PlotDotPlots.R"))
source(file.path(base_dir, "scripts", "PlotBoxPlots.R"))
source(file.path(base_dir, "scripts", "PlotFacetBoxPlots.R"))


# Analysis options ------------------------------------------------------------

save_csv               <- FALSE # Set to TRUE to save descriptive statistics tables as CSV files in the derivatives directory
show_pie_charts        <- FALSE # Set to TRUE to create pie charts of categorical variables (e.g., species
show_amcharts          <- FALSE # Set to TRUE to create interactive 3D pie charts using amCharts4
show_dot_plots         <- TRUE # Set to TRUE to create dot plots of spectral quality metrics by different grouping variables
show_box_plots         <- TRUE # Set to TRUE to create box plots of spectral quality metrics by different grouping variables
show_facet_plots       <- TRUE # Set to TRUE to create facet plots of spectral quality metrics by different grouping variables
show_model_diagnostics <- TRUE # Set to TRUE to show model diagnostic plots (e.g., residuals, Q-Q plots) for linear mixed-effects models
calc_VPCs              <- TRUE # Set to TRUE to calculate variance partition coefficients (VPCs) from linear mixed-effects models to assess the proportion of variance explained by each random effect
export_model_table     <- TRUE # Set to TRUE to export a modelsummary comparison table of LMEM results to Word (.docx)
run_LRTs               <- TRUE # Set to TRUE to run likelihood ratio tests (LRTs) to compare linear mixed-effects models with different random effects structures and derive p-values for the added random effects
run_pbkrtest           <- FALSE # Set to TRUE to run parametric bootstrapping using the pbkrtest package to compare linear mixed-effects models with different random effects structures and derive p-values for the added random effects (can be time-consuming with larger datasets)

# Load data -------------------------------------------------------------------
# Also clean up data (incl. outlier removal) and create new variables (e.g., normalized SNR/LW ratio)

DATA <- LoadData(
  csv_file = file.path(data_dir, "CoMP_MRS_Rstats_input.csv"),
  outl_rm_strategy = "group",
  verbose = TRUE
)


# Run descriptive statistics --------------------------------------------------

STATS <- DescripStats(data = DATA$data)


# Save descriptive statistics tables as CSV files -----------------------------

if (save_csv) {
  SaveCSV(data = STATS, out_dir = csv_dir)
}


# Plot pie charts -------------------------------------------------------------

if (show_pie_charts) {
  PlotPieCharts(
    data = DATA$data,
    out_dir = file.path(plots_dir, "pie_charts"),
    plots_dir = plots_dir,
    show_pie_charts = show_pie_charts,
    show_amcharts = show_amcharts
  )
}


# Plot dot plots --------------------------------------------------------------

if (show_dot_plots) {
  PlotDotPlots(
    data = DATA$data,
    data_orig = DATA$data_orig,
    out_dir = file.path(plots_dir, "dot_plots")
  )
}


# Plot box plots --------------------------------------------------------------

if (show_box_plots) {
  y_vars <- list(
    #list(var = "LW_norm",             label = "Normalized LW"),
    #list(var = "SNR_norm",            label = "Normalized SNR"),
    #list(var = "SNR_LW_Product_norm", label = "Normalized SNR×LW product"),
    list(var = "SNR_LW_Ratio_norm", label = "Normalized SNR/LW ratio")
  )
  
  x_vars <- c(
    "DP",
    "SiteID",
    "Species",
    "Sex",
    "Vendor",
    "FieldStrength",
    "Sequence",
    "Sequence_collapsed",
    "VOI",
    "Cryoprobe",
    "ShimMethod"
  )
  
  DATA$data <- DATA$data %>%
    dplyr::mutate(
      VOI = stringr::str_replace_all(
        VOI,
        c(
          "Rhippocampus" = "R-hipp",
          "Lhippocampus" = "L-hipp",
          "Rstriatum"    = "R-str",
          "Lstriatum"    = "L-str"
        )
      )
    )
  
  all_boxplots <- PlotBoxPlots(
    data = DATA$data,
    out_dir = file.path(plots_dir, "box_plots"),
    y_vars  = y_vars,
    x_vars  = x_vars
  )
}


# Plot faceted box plots ------------------------------------------------------

if (show_facet_plots) {
  y_vars <- list(
    # list(var = "LW_norm",             label = "Normalized LW"),
    # list(var = "SNR_norm",            label = "Normalized SNR"),
    # list(var = "SNR_LW_Product_norm", label = "Normalized SNR×LW product"),
    list(var = "SNR_LW_Ratio_norm",   label = "Normalized SNR/LW ratio")
  )
  
  x_vars <- list(
    list(var = "DP", label = "Data Packet"),
    list(var = "SiteID", label = "Site ID"),
    list(var = "Species", label = "Species"),
    list(var = "Cryoprobe", label = "Cryoprobe"),
    list(var = "ShimMethod", label = "Shim method"),
    list(var = "FieldStrength", label = "Field strength")
  )
  
  facet_vars <- list(
    list(var = "Sequence", label = "MRS sequence")
  )
  
  facet_plots <- PlotFacetBoxPlots(
    data      = DATA$data,
    out_dir   = file.path(plots_dir, "facet_box_plots"),
    x_var     = x_vars,
    y_var     = y_vars,
    facet_var = facet_vars
  )
}


# Run linear mixed-effects modeling -------------------------------------------

source(file.path(base_dir, "scripts", "RunLMEM.R"))

LMEM_MODELS <- list() # Initialize list to store LMEM models

# Variables

dv <- "SNR_LW_Ratio_norm"

random_effects <- list(
  M.0 = list(
    DP = "1",
    SiteID = "1",
    Vendor = "1"
  ),
  M.0.a = list(
    Vendor = "1"
  ),
  M.0.b = list(
    SiteID = "1"
  ),
  M.0.c = list(
    DP = "1"
  ),
  M.0.d = list(
    DP = "1",
    SiteID = "1"
  ),
  M.1 = list(
    DP = "1",
    SiteID = "1",
    Vendor = "1",
    Species = "1"
  ),
  M.1.a = list(
    DP = "1",
    ShimMethod = "1"
  ),
  M.1.b = list(
    DP = "1",
    Cryoprobe = "1"
  ),
  M.1.c = list(
    DP = "1",
    Species = "1"
  ),
  M.1.d = list(
    DP = "1",
    VOI = "1"
  ),
  M.1.e = list(
    DP = "1",
    Sequence = "1"
  ),
  M.1.f = list(
    DP = "1"
  ),
  M.1.g = list(
    DP = "1"
  ),
  M.1.h = list(
    DP = "1"
  ),
  M.2 = list(
    DP = "1",
    Sequence = "1",
    Cryoprobe = "1"
  )
)

fixed_effects <- list(
  M.1.f = c("FieldStrength"),
  M.1.g = c("Age"),
  M.1.h = c("Sex")
)

# # Null model with no random effects
# 
# M.null <- lm(
#   formula = as.formula(paste(dv, "~ 1")),
#   data = DATA$data
# )

# Run LMEM models

LMEM_MODELS <- lapply(names(random_effects), function(model_name) {
  
  if (is_empty(fixed_effects[[model_name]])) {
    fix_ef <- ""
  } else {
    fix_ef = fixed_effects[[model_name]]
  }
  
  RunLMEM(
    data = DATA$data,
    dv = dv,
    rand_ef = random_effects[[model_name]],
    fix_ef = fix_ef
  )
  
})

names(LMEM_MODELS) <- names(random_effects)

### Model diagnostics ---------------------------------------------------------

if (show_model_diagnostics) {
  
  LMEM_MODELS.diagnostics <- lapply(
    names(random_effects),
    function(model_name) {
      model <- LMEM_MODELS[[model_name]]
      diagnostics <- list(
        random_effects      = ranef(model, condVar = TRUE),
        random_effects.plot = dotplot(ranef(model, condVar = TRUE), scales = "free"),
        model_performance   = model_performance(model),
        check_model         = check_model(model),
        check_singularity   = check_singularity(model)
      )
      return(diagnostics)
    }
  )
  
  names(LMEM_MODELS.diagnostics) <- names(random_effects)
  
}

### Variance partitioning -----------------------------------------------------

source(file.path(base_dir, "scripts", "ExtractVPCs.R"))

if (calc_VPCs) {
  VPCs <- ExtractVPCs(LMEM_MODELS, verbose = FALSE)
}

### Inference by LRT with parametric bootstrapping ----------------------------
# Compare large model with smaller model to derive p-value for added random 
# effect (e.g., Sequence)
# Note, models have to be fitted with REML = FALSE for valid comparison by LRT,
# and this can be time-consuming with larger datasets

source(file.path(base_dir, "scripts", "RunLRT.R"))

model_contrasts <- list(
  small_models = c("M.0.b", "M.0.c", "M.0.c", "M.0.c", "M.1.e"),
  large_models = c("M.0.d", "M.0.d", "M.1.b", "M.1.e", "M.2")
)

if (run_LRTs || run_pbkrtest) {
  LRT_RESULTS <- RunLRT(
    models          = LMEM_MODELS,
    model_contrasts = model_contrasts,
    run_LRTs        = run_LRTs,
    run_pbkrtest    = run_pbkrtest,
    verbose         = FALSE
  )
}


### Export LMEM model comparison table to Word --------------------------------

source(file.path(base_dir, "scripts", "ExportModelTable_old.R"))
source(file.path(base_dir, "scripts", "ExportModelTable.R"))

if (export_model_table) {
  # MODEL_TABLE.old <- ExportModelTable_old(
  #   # models = c(list(M.null = M.null), LMEM_MODELS),
  #   models = LMEM_MODELS,
  #   vpcs = if (calc_VPCs) VPCs else NULL,
  #   out_dir = deriv_dir
  # )
  
  MODEL_TABLE <- ExportModelTable(
    # models = c(list(M.null = M.null), LMEM_MODELS),
    models = LMEM_MODELS,
    vpcs = if (calc_VPCs) VPCs else NULL,
    out_dir = deriv_dir
  )
}
