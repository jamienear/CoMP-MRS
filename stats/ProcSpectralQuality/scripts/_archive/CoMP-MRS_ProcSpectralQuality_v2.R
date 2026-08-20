# R script to visualize and statistically analyze data from CoMP-MRS
#
# Authors:
#   Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
#   Diana G. Rotaru, Ph.D. (diana.rotaru@meduniwien.ac.at)
# Last updated: 2026-03-22


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

cv.bs <- function(x) {
  return(apply(x, 2, function(y) sd(y, na.rm = TRUE)) / colMeans(x, na.rm = TRUE) * 100)
}

cv.bs.all <- function(x) {
  return(sd(as.numeric(x), na.rm = TRUE) / mean(as.numeric(x), na.rm = TRUE) * 100)
}

cv.ws <- function(x) {
  return(apply(x, 1, function(y) sd(y, na.rm = TRUE)) / rowMeans(x, na.rm = TRUE) * 100)
}

cv.ws.dpm <- function(x) {
  # Calculate the within-subject CV using Åsberg and Bolann's difference in
  # percentage of the mean approach (Åsberg & Bolann. doi:10.1093/jalm/jfac098).
  d <- x[, 2] - x[, 1]
  m <- (x[, 1] + x[, 2]) / 2
  dpm <- (d / m) * 100
  sd_dpm <- sd(dpm, na.rm = TRUE)
  return(sd_dpm / sqrt(2))
}

# This function prints all of the descriptive stats in a list
print.descrip.stats <- function(x) {
  list(
    mean.by_ses       = round(colMeans(x, na.rm = TRUE), 2),
    mean.overall      = round(mean(as.numeric(x), na.rm = TRUE), 2),
    std.by_ses        = round(apply(x, 2, function(x) sd(x, na.rm = TRUE)), 2),
    std.overall       = round(sd(as.numeric(x), na.rm = TRUE), 2),
    cv.bs.by_ses      = round(cv.bs(x), 1),
    cv.ws.by_subj     = round(cv.ws(x), 1),
    # cv.ws.overall     = round(mean(cv.ws(x), na.rm = TRUE), 1),
    cv.ws.dpm.overall = round(cv.ws.dpm(x), 1),
    cv.bs.overall     = round(cv.bs.all(x), 1),
    ICC.psych         = ICC(x),
    ICC.irr           = icc(x, model = "twoway", type = "agreement", unit = "single")
  )
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

# R-squared functions for linear mixed models (Nakagawa & Schielzeth, 2013, doi:10.1111/j.2041-210X.2012.00261.x)
r2.global <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

r2.local <- function(small_model, large_model) {
  vc_small_model <- as.data.frame(VarCorr(small_model))[4]
  vc_large_model <- as.data.frame(VarCorr(large_model))[4]
  out <- as.numeric((tail(vc_small_model, n=1) - tail(vc_large_model, n=1)) / tail(vc_small_model, n=1))
  return(out)
}


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
    AnimalAge = col_double(),
    AnimalSex = col_factor(),
    AnimalWeight = col_double(),
    MRvendor = col_factor(),
    MRfield = col_double(),
    MRsequence = col_factor(),
    MRbrainregion = col_factor(),
    MRvoxelvolume = col_double(),
    MRaverages = col_double(),
    LW = col_double(),
    SNR = col_double(),
    SNR_LW_Ratio = col_double(),
    CompCheck = col_factor()
  )
)

DATA[is.na(DATA)] <- NA


# # Pad data frames with NA's when there are missing datasets -------------------
# 
# empty_rows <- list(
#   "hermes_press"    = data.frame(ind = c(27, 28)),
#   "hermes_slaser"   = data.frame(ind = c(27, 28)),
#   "hercules_press"  = data.frame(ind = c(24)),
#   "hercules_slaser" = data.frame(ind = c(16, 25, 26))
# )
# 
# for (ii in acq) {
#   DATA_diff1[[ii]] <- insertRows(DATA_diff1[[ii]], r = empty_rows[[ii]]$ind, rcurrent = F)
#   DATA_diff2[[ii]] <- insertRows(DATA_diff2[[ii]], r = empty_rows[[ii]]$ind, rcurrent = F)
#   DATA_sum[[ii]] <- insertRows(DATA_sum[[ii]], r = empty_rows[[ii]]$ind, rcurrent = F)
# 
#   DATA_diff1[[ii]]$sub <- as_factor(subs)
#   DATA_diff1[[ii]]$ses <- as_factor(sess)
#   DATA_diff1[[ii]]$acq <- as_factor(ii)
#   DATA_diff2[[ii]]$sub <- as_factor(subs)
#   DATA_diff2[[ii]]$ses <- as_factor(sess)
#   DATA_diff2[[ii]]$acq <- as_factor(ii)
#   DATA_sum[[ii]]$sub <- as_factor(subs)
#   DATA_sum[[ii]]$ses <- as_factor(sess)
#   DATA_sum[[ii]]$acq <- as_factor(ii)
# }


# Data normalization -----------------------------------------------------------

DATA$SNR_LW_Ratio_norm <- 
  DATA$SNR_LW_Ratio/(sqrt(DATA$MRaverages)*DATA$MRvoxelvolume)


# Descriptive stats -----------------------------------------------------------

STATS <- list(
  mean = list(),
  sd = list(),
  cv = list()
)

STATS$mean$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(mean = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$mean$AnimalSpecies$Uncorr <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  summarise(mean = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$mean$AnimalSpecies$Norm <-
  DATA %>%
  group_by(MRvendor, AnimalSpecies) %>%
  #  mutate(SNR_LW_Ratio_norm = SNR_LW_Ratio/(sqrt(MRaverages)*MRvoxelvolume)) %>%
  summarise(mean = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$mean$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(mean = mean(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sd$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(sd = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$sd$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(sd = sd(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cv$MRvendor <-
  DATA %>%
  group_by(MRvendor) %>%
  summarise(cv = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cv$DP <-
  DATA %>%
  group_by(DP) %>%
  summarise(cv = cv(SNR_LW_Ratio_norm, na.rm = TRUE))

STATS$cv$subj <- cv(DATA$SNR_LW_Ratio_norm, na.rm = TRUE)

# Quality metrics
# quality.metrics <- summary(DATA)


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


# Linear-mixed models ---------------------------------------------------------

# Selection of optimizers to use
defaultControl <- list(algorithm="NLOPT_LN_BOBYQA", xtol_rel=1e-6, maxeval=1e5)
optWrap.nloptwrap.NLOPT_LN_BOBYQA <- lmerControl(optimizer="nloptFun")
optWrap.optimx.nlminb             <- lmerControl(optimizer="optimx", optCtrl=list(method="nlminb", eval.max=1e5))
optWrap.minqa.bobyqa              <- lmerControl(optimizer="bobyqa")

optimToUse <- optWrap.minqa.bobyqa

Y.M.0.0  <- lm(SNR_LW_Ratio_norm ~ 1, data=DATA)
Y.M.0.1 <- lmer(SNR_LW_Ratio_norm ~ (1 | DP), data=DATA, REML=F, control=optimToUse)
Y.M.0.2 <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor), data=DATA, REML=F, control=optimToUse)
Y.M.0.3 <- lmer(SNR_LW_Ratio_norm ~ (1 | SiteID), data=DATA, REML=F, control=optimToUse)
Y.M.0.4  <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) , data=DATA, REML=F, control=optimToUse)
Y.M.0.5  <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) , data=DATA, REML=F, control=optimToUse)
# Y.M.1.0  <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) , data=DATA, REML=F, control=optimToUse)
# Y.M.2.0  <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) , data=DATA, REML=F, control=optimToUse)
# Y.M.3.0  <- lmer(SNR_LW_Ratio_norm ~ (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data=DATA, REML=F, control=optimToUse)
# Y.M.4.0  <- lmer(SNR_LW_Ratio_norm ~ MRfield + (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data=DATA, REML=F, control=optimToUse)
# Y.M.5.0  <- lmer(SNR_LW_Ratio_norm ~ MRfield + AnimalAge + (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data=DATA, REML=F, control=optimToUse)
# Y.M.6.0  <- lmer(SNR_LW_Ratio_norm ~ MRfield + AnimalAge + AnimalSex + (1 | MRvendor) + (1 | SiteID) + (1 | DP) + 
#                    (1 | AnimalSpecies) + (1 | MRbrainregion) + (1 | MRsequence) , data=DATA, REML=F, control=optimToUse)

# Model diagnostics ------------------------------------------------------

# check_model(Y.M.0.0)
check_model(Y.M.0.1)
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


# Variance components -----------------------------------------------------

# Save variance components of three-level null model
print(VarCorr(Y.M.0.1), comp = c("Variance", "Std.Dev."))
vc <- as.data.frame(VarCorr(Y.M.0.1))[4]
vpc.M.0.1 <- round(vc/sum(vc)*100,1)
dimnames(vpc.M.0.1)[[1]] <- c("DP", "Residual")



# NAAG_tCr vs. Asp_tCr --------------------------------------------------------

data.merged <- rbind(
  DATA_sum$hercules_press,
  DATA_sum$hercules_slaser
)

# MMCD outlier rejection
outl <- MMCD(data.merged$NAAG_tCr, data.merged$Asp_tCr)

data.merged.no_outl <- data.merged[outl, ]

ggplot(data.merged.no_outl, aes(x = NAAG_tCr, y = Asp_tCr)) +
  geom_point() +
  geom_smooth(method = lm, color = "red", fill = "gray", se = TRUE) +
  theme_classic2()

lm.out <- lm(NAAG_tCr ~ Asp_tCr, data.merged.no_outl)
summary(lm.out)
