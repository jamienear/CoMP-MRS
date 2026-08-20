# R script to visualize and statistically analyze data from CoMP-MRS
#
# Authors:
#   Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
#   Diana G. Rotaru, Ph.D. (diana.rotaru@meduniwien.ac.at)
#
# Last updated: 2026-03-23


# Initialize ------------------------------------------------------------------

rm(list = ls()) # Clear the current environment
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE) # If using RStudio, this clears all plots
cat("\014") # CTRL+L (clear console)


# Set data directories --------------------------------------------------------

##### CHANGE THESE DIRECTORIES AS NEEDED #####

# Set data directories
curr_dir  <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
data_dir  <- file.path(curr_dir, "data")
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
  "performance"
)

load.packages(packages)


# Set inline functions --------------------------------------------------------

# MAD outlier (Leys et al., 2013, doi:10.1016/j.jesp.2013.03.013)
mad.outlier <- function(x, thresh = 2.5) {
  m <- mad(x, na.rm = TRUE)
  lb <- median(x, na.rm = TRUE) - thresh * m
  ub <- median(x, na.rm = TRUE) + thresh * m
  outl <- x < lb | x > ub
  return(outl)
}

# Removal of multivariate outliers using the Mahalanobis-minimum covariance 
# determinant (MMCD) (Leys et al., 2018, doi:10.1016/j.jesp.2017.09.011)
MMCD <- function(x, y, quant = 0.75, alpha = 0.01) {
  mat <- cbind(x, y)
  output <- cov.mcd(na.omit(mat), quantile.used = nrow(na.omit(mat)) * quant, nsamp = "exact")
  mmcd <- mahalanobis(mat, output$center, output$cov)
  cutoff <- qchisq(p = 1 - alpha, df = 2)
  out <- -which(mmcd > cutoff | is.na(mmcd))
  if (length(out) == 0) {
    out <- 1:nrow(mat)
  }
  return(out)
}

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

# Load the participants.tsv file
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
    MRcoil = col_factor(),
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

DATA[is.na(DATA)] <- NA # Ensure that all missing values are represented as NA
                        # (in case there are any other representations of missing
                        # data in the original CSV)



# Spectral quality metrics and normalization ----------------------------------------------------------
# Get both the ratio and the product of SNR and LW and the normalize 
#   the values to the voxel volume and square root of avereges number

DATA$SNR_LW_Ratio_norm <- 
  DATA$SNR_LW_Ratio / (sqrt(DATA$MRaverages) * DATA$MRvoxelvolume)

DATA$SNR_LW_Product_norm <-
  (DATA$LW * DATA$SNR)/ (sqrt(DATA$MRaverages) * DATA$MRvoxelvolume)

DATA$LW_norm <- 
  DATA$LW / (sqrt(DATA$MRaverages) * DATA$MRvoxelvolume)

DATA$SNR_norm <- 
  DATA$SNR / (sqrt(DATA$MRaverages) * DATA$MRvoxelvolume)

# Descriptive stats -----------------------------------------------------------

STATS <- list(
  mean = list(),
  sd = list(),
  cv = list()
)

### DP ---------------------------------------------------------

STATS$meanLW$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$cvLW$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$meanSNR$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$cvSNR$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$meanRatio$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvRatio$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvProduct$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### SiteID ---------------------------------------------------------

STATS$meanLW$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$SiteID <-
  DATA %>%
  group_by(SiteID) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### AnimalSpecies ---------------------------------------------------------

STATS$meanLW$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$AnimalSpecies <-
  DATA %>%
  group_by(AnimalSpecies) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### MRvendor ---------------------------------------------------------

STATS$meanLW$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### MRfield ---------------------------------------------------------

STATS$meanLW$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$MRfield <-
  DATA %>%
  group_by(MRfield) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### MRsequence ---------------------------------------------------------

STATS$meanLW$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$MRsequence <-
  DATA %>%
  group_by(MRsequence) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### MRbrainregion ---------------------------------------------------------

STATS$meanLW$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$MRbrainregion <-
  DATA %>%
  group_by(MRbrainregion) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### AnimalSpecies x MRvendor --------------------------------------------------

STATS$meanLW$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$AnimalSpeciesMRvendor <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

### AnimalSpecies x MRbrainregion ---------------------------------------------

STATS$meanLW$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(meanLW = mean(LW_norm, na.rm = TRUE))

STATS$sdLW$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(sdLW = sd(LW_norm, na.rm = TRUE))

STATS$meanSNR$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(meanSNR = mean(SNR_norm, na.rm = TRUE))

STATS$sdSNR$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(sdSNR = sd(SNR_norm, na.rm = TRUE))

STATS$meanRatio$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(meanRatio = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sdRatio$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(sdRatio = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$meanProduct$AnimalSpeciesMRbrainregion<-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(meanProduct = mean(SNR_LW_Product_norm, na.rm = TRUE))

STATS$sdProduct$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(sdProduct = sd(SNR_LW_Product_norm, na.rm = TRUE))

STATS$cvLW$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(cvLW = cv(LW_norm, na.rm = TRUE))

STATS$cvSNR$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(cvSNR = cv(SNR_norm, na.rm = TRUE))

STATS$cvRatio$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(cvRatio = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cvProduct$AnimalSpeciesMRbrainregion <-
  DATA %>%
  group_by(AnimalSpecies, MRbrainregion) %>%
  summarise(cvProduct = cv(SNR_LW_Product_norm, na.rm = TRUE))

## If data not already normalized, use mutate
# STATS$mean$AnimalSpecies$Norm <-
# DATA %>%
# group_by(MRvendor, AnimalSpecies) %>%
# mutate(SNR_LW_Ratio_norm = SNR_LW_Ratio/(sqrt(MRaverages)*MRvoxelvolume)) %>%
# summarise(mean = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

# STATS$mean$AnimalSpecies$Norm <-
  # DATA %>%
  # group_by(MRvendor, AnimalSpecies) %>%
  #  mutate(SNR_LW_Ratio_norm = SNR_LW_Ratio/(sqrt(MRaverages)*MRvoxelvolume)) %>%
  # summarise(mean = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

### All -----------------------------------------------------------------------

STATS$mean$allLW <- mean(DATA$LW_norm, na.rm = TRUE)
STATS$sd$allLW <- sd(DATA$LW_norm, na.rm = TRUE)
STATS$cv$allLW <- cv(DATA$LW_norm, na.rm = TRUE)

STATS$mean$allSNR <- mean(DATA$SNR_norm, na.rm = TRUE)
STATS$sd$allSNR <- sd(DATA$SNR_norm, na.rm = TRUE)
STATS$cv$allSNR <- cv(DATA$SNR_norm, na.rm = TRUE)

STATS$mean$allRatio <- mean(DATA$SNR_LW_Ratio_norm, na.rm = TRUE)
STATS$sd$allRatio <- sd(DATA$SNR_LW_Ratio_norm, na.rm = TRUE)
STATS$cv$allRatio <- cv(DATA$SNR_LW_Ratio_norm, na.rm = TRUE)

STATS$mean$allProduct <- mean(DATA$SNR_LW_Product_norm, na.rm = TRUE)
STATS$sd$allProduct <- sd(DATA$SNR_LW_Product_norm, na.rm = TRUE)
STATS$cv$allProduct <- cv(DATA$SNR_LW_Product_norm, na.rm = TRUE)

#  Pie charts -----------------------------------------------------------------

ggplot(DATA, aes(x="", y=AnimalSpecies, fill=AnimalSpecies)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
    plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
    legend.position = "bottom",
    legend.background = element_blank(),
    legend.title = element_text(face = "bold", size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    text = element_text(family = "Helvetica"),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", lineend = "square", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25, lineend = "square")
  )

ggsave(
  file.path(plots_dir, "AnimalSpiecies_piechart.png"),
  height = 5, width = 5, units = "in"
)


# Boxplots --------------------------------------------------------------------

bwplot(
  SNR_LW_Ratio_norm ~ MRvendor,
  data = DATA,
  xlab = "Vendor",
  ylab = "SNR_LW_Ratio_norm"
)

bwplot(
  SNR_LW_Ratio_norm ~ MRvendor | AnimalSpecies,
  data = DATA,
  xlab = "Vendor",
  ylab = "SNR_LW_Ratio_norm"
)

bwplot(
  SNR_LW_Ratio_norm ~ SiteID,
  data = DATA,
  xlab = "Site",
  ylab = "SNR_LW_Ratio_norm"
)

bwplot(
  SNR_LW_Ratio_norm ~ DP,
  data = DATA,
  xlab = "DP",
  ylab = "SNR_LW_Ratio_norm"
)


# Linear-mixed models ---------------------------------------------------------

# Selection of optimizers to use
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA", xtol_rel=1e-6, maxeval=1e5)
optWrap.nloptwrap.NLOPT_LN_BOBYQA <- lmerControl(optimizer="nloptFun")
optWrap.optimx.nlminb             <- lmerControl(optimizer="optimx", optCtrl=list(method="nlminb", eval.max=1e5))
optWrap.minqa.bobyqa              <- lmerControl(optimizer="bobyqa")

optimToUse <- optWrap.nloptwrap.NLOPT_LN_BOBYQA

# Y.M.0.0 <- lm(SNR_LW_Ratio_norm ~ 1, data=DATA)
Y.M.0.1 <- lmer(cs.(SNR_LW_Ratio_norm) ~ (1 | DP), data = DATA, REML = FALSE, control = optimToUse)
# Y.M.0.2 <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor), data = DATA, REML = FALSE, control = optimToUse)
# Y.M.0.3 <- lmer(SNR_LW_Ratio_norm ~ (1 | SiteID), data = DATA, REML = FALSE, control = optimToUse)
# Y.M.0.4 <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) , data = DATA, REML = FALSE, control = optimToUse)
Y.M.0.5 <- lmer(cs.(SNR_LW_Ratio_norm) ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP),
                data = DATA, REML = FALSE, control = optimToUse)
# Y.M.1.0 <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) , data = DATA, REML = FALSE, control = optimToUse)
# Y.M.2.0 <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) , data = DATA, REML = FALSE, control = optimToUse)
# Y.M.3.0 <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data = DATA, REML = FALSE, control = optimToUse)
# Y.M.4.0 <- lmer(SNR_LW_Ratio_norm ~ MRfield + (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data = DATA, REML = FALSE, control = optimToUse)
# Y.M.5.0 <- lmer(SNR_LW_Ratio_norm ~ MRfield + AnimalAge + (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data = DATA, REML = FALSE, control = optimToUse)
# Y.M.6.0 <- lmer(SNR_LW_Ratio_norm ~ MRfield + AnimalAge + AnimalSex + (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data = DATA, REML = FALSE, control = optimToUse)


# Model diagnostics -----------------------------------------------------------

dotplot(ranef(Y.M.0.1, condVar=TRUE), scales=list(relation="free"))
dotplot(ranef(Y.M.0.5, condVar=TRUE), scales=list(relation="free"))

# check_model(Y.M.0.0)
check_model(Y.M.0.1)
check_model(Y.M.0.5)
model_performance(Y.M.0.1)
model_performance(Y.M.0.5)
# r2(Y.M.0.1)
# check_model(Y.M.0.2)
# check_model(Y.M.0.3)
# check_model(Y.M.0.4)
# check_model(Y.M.0.5)
# check_model(Y.M.1.0)
# check_model(Y.M.2.0)
# check_model(Y.M.3.0)
# check_model(Y.M.4.0)
# check_model(Y.M.5.0)
# check_model(Y.M.6.0)

# z <- profile(Y.M.0.1)
# z_ci <- confint(z, level=0.95)^2
# xyplot(z, aspect=1)

# plot(Y.M.0.1)
# hist(resid(Y.M.0.1))


# Variance components ---------------------------------------------------------

# Save variance components of three-level null model
# print(VarCorr(Y.M.0.1), comp = c("Variance", "Std.Dev."))
vc <- as.data.frame(VarCorr(Y.M.0.1))[4]
vpc.M.0.1 <- round(vc/sum(vc)*100,1)
dimnames(vpc.M.0.1)[[1]] <- c("DP", "Residual")

# print(VarCorr(Y.M.0.5), comp = c("Variance", "Std.Dev."))
vc <- as.data.frame(VarCorr(Y.M.0.5))[4]
vpc.M.0.5 <- round(vc/sum(vc)*100,1)
dimnames(vpc.M.0.5)[[1]] <- c("DP", "Site", "Vendor", "Residual")



