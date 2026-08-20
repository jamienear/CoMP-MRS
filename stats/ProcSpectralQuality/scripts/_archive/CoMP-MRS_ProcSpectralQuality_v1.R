# R script to visualize and statistically analyze data from "Test-retest
# reliability of multi-metabolite edited MRS at 3T using PRESS and sLASER"
#
# Author: Mark Mikkelsen, Ph.D. (mam4041@med.cornell.edu)
# Last updated: 2025-10-22


# Initialize ------------------------------------------------------------------

rm(list = ls()) # Clear the current environment
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE) # If using RStudio, this clears all plots
cat("\014") # CTRL+L (clear console)


# Set data directories --------------------------------------------------------

##### CHANGE THESE DIRECTORIES AS NEEDED #####

# Set data directories
curr_dir <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
data_dir <- file.path(curr_dir, "data")
deriv_dir <- file.path(data_dir, "derivatives")
out_dir <- file.path(deriv_dir, "output")
BA_plots <- file.path(out_dir, "BA_plots")

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

if (!dir.exists(BA_plots)) {
  dir.create(BA_plots)
} else {
  unlink(BA_plots, recursive = TRUE)
  dir.create(BA_plots)
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
  "latex2exp"
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


# Load data -------------------------------------------------------------------

# Load the participants.tsv file
data <- read_csv(
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

data[is.na(data)] <- NA



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


# Descriptive stats -----------------------------------------------------------

qualitymetrics <- summary(data)

#  Pie charts -----------------------------------------------------------------

ggplot(data, aes(x="", y=AnimalSpecies, fill=AnimalSpecies)) +
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
    text = element_text(family = "Arial"),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", lineend = "square", linewidth = 0.25),
    axis.ticks = element_line(color = "black", linewidth = 0.25, lineend = "square")
  ) 
ggsave(
  file.path(out_dir, "AnimalSpecies_piechart.png"),
  height = 5, width = 5, units = "in"
)

# Bland-Altman plots ----------------------------------------------------------

xLims <- list(
  "Asc_tCr" = data.frame(
    hermes = c(0, 0.12, 0.02),
    hercules = c(0, 0.3, 0.05)
  ),
  "Asp_tCr" = data.frame(
    hermes = c(0, 0.25, 0.05),
    hercules = c(0.15, 0.5, 0.05)
  ),
  "GABA_tCr" = data.frame(
    hermes = c(0.15, 0.35, 0.05),
    hercules = c(0.15, 0.275, 0.05)
  ),
  "Gln_tCr" = data.frame(
    hermes = c(0, 0.5, 0.1),
    hercules = c(0, 0.5, 0.1)
  ),
  "Glu_tCr" = data.frame(
    hermes = c(0.5, 1.75, 0.25),
    hercules = c(0.5, 1.75, 0.25)
  ),
  "GSH_tCr" = data.frame(
    hermes = c(0.15, 0.325, 0.05),
    hercules = c(0.15, 0.325, 0.05)
  ),
  "Lac_tCr" = data.frame(
    hermes = c(0, 0.5, 0.1),
    hercules = c(0.025, 0.175, 0.05)
  ),
  "NAA_tCr" = data.frame(
    hermes = c(1, 1.75, 0.25),
    hercules = c(1, 1.75, 0.25)
  ),
  "NAAG_tCr" = data.frame(
    hermes = c(0, 0.2, 0.05),
    hercules = c(0.05, 0.25, 0.05)
  ),
  "mI_tCr" = data.frame(
    hermes = c(0.5, 1.25, 0.25),
    hercules = c(0.5, 1.25, 0.25)
  ),
  "tCho_tCr" = data.frame(
    hermes = c(0.15, 0.25, 0.025),
    hercules = c(0.15, 0.25, 0.025)
  )
)

yLims <- list(
  "Asc_tCr" = data.frame(
    hermes = c(-450, 450, 150),
    hercules = c(-250, 250, 50)
  ),
  "Asp_tCr" = data.frame(
    hermes = c(-450, 450, 150),
    hercules = c(-150, 150, 50)
  ),
  "GABA_tCr" = data.frame(
    hermes = c(-100, 100, 25),
    hercules = c(-75, 75, 25)
  ),
  "Gln_tCr" = data.frame(
    hermes = c(-150, 150, 50),
    hercules = c(-150, 150, 50)
  ),
  "Glu_tCr" = data.frame(
    hermes = c(-30, 30, 10),
    hercules = c(-50, 50, 10)
  ),
  "GSH_tCr" = data.frame(
    hermes = c(-100, 100, 20),
    hercules = c(-150, 150, 50)
  ),
  "Lac_tCr" = data.frame(
    hermes = c(-350, 350, 50),
    hercules = c(-100, 100, 25)
  ),
  "NAA_tCr" = data.frame(
    hermes = c(-50, 50, 10),
    hercules = c(-20, 20, 5)
  ),
  "NAAG_tCr" = data.frame(
    hermes = c(-300, 300, 50),
    hercules = c(-300, 300, 100)
  ),
  "mI_tCr" = data.frame(
    hermes = c(-25, 25, 5),
    hercules = c(-40, 40, 10)
  ),
  "tCho_tCr" = data.frame(
    hermes = c(-25, 25, 5),
    hercules = c(-50, 50, 10)
  )
)

for (ii in subspec) {
  for (kk in metab_tCr) {
    if (
      ((ii == "DATA_diff2") && (kk == "Asc_tCr")) ||
        ((ii == "DATA_diff2") && (kk == "Asp_tCr")) ||
        ((ii == "DATA_diff1") && (kk == "GABA_tCr")) ||
        ((ii == "DATA_sum") && (kk == "Gln_tCr")) ||
        ((ii == "DATA_sum") && (kk == "Glu_tCr")) ||
        ((ii == "DATA_diff2") && (kk == "GSH_tCr")) ||
        ((ii == "DATA_diff2") && (kk == "Lac_tCr")) ||
        ((ii == "DATA_diff2") && (kk == "NAAG_tCr")) ||
        ((ii == "DATA_sum") && (kk == "NAA_tCr")) ||
        ((ii == "DATA_sum") && (kk == "mI_tCr")) ||
        ((ii == "DATA_sum") && (kk == "tCho_tCr"))
    ) {
      for (jj in 1:2) {
        # Skip these metabolites for HERMES data
        if (any(kk == c("Asc_tCr", "Asp_tCr", "Lac_tCr", "NAAG_tCr")) && jj == 1) {
          next
        }

        tmp_press_ses01 <- subset(eval(parse(text = paste0(ii, "$", scheme[jj], "_press"))), 
                                  acq == paste0(scheme[jj], "_press") & ses == "ses-01")
        tmp_press_ses02 <- subset(eval(parse(text = paste0(ii, "$", scheme[jj], "_press"))), 
                                  acq == paste0(scheme[jj], "_press") & ses == "ses-02")
        tmp_slaser_ses01 <- subset(eval(parse(text = paste0(ii, "$", scheme[jj], "_slaser"))), 
                                   acq == paste0(scheme[jj], "_slaser") & ses == "ses-01")
        tmp_slaser_ses02 <- subset(eval(parse(text = paste0(ii, "$", scheme[jj], "_slaser"))), 
                                   acq == paste0(scheme[jj], "_slaser") & ses == "ses-02")

        # Because Asc was not detectable in HERCULES-sLASER data and so the 
        # values are all 0's, create an array of random numbers to prevent the 
        # MMCD outlier rejection subroutine from failing
        if (ii == "DATA_diff2" && kk == "Asc_tCr" && jj == 2) {
          tmp_slaser_ses01$Asc_tCr <- rnorm(nrow(tmp_slaser_ses01),
            mean = mean(tmp_press_ses01$Asc_tCr, na.rm = TRUE),
            sd = sd(tmp_press_ses01$Asc_tCr, na.rm = TRUE)
          )
          tmp_slaser_ses02$Asc_tCr <- rnorm(nrow(tmp_slaser_ses02),
            mean = mean(tmp_press_ses02$Asc_tCr, na.rm = TRUE),
            sd = sd(tmp_press_ses02$Asc_tCr, na.rm = TRUE)
          )
        }

        tmp <- as_tibble(data.frame(
          sub    = as.factor(c(tmp_press_ses01$sub, tmp_slaser_ses01$sub)),
          acqs   = as.factor(c(rep(acq[jj * 2 - 1], times = nrow(tmp_press_ses01)), rep(acq[jj * 2], times = nrow(tmp_slaser_ses01)))),
          ses_01 = c(tmp_press_ses01[[kk]], tmp_slaser_ses01[[kk]]),
          ses_02 = c(tmp_press_ses02[[kk]], tmp_slaser_ses02[[kk]])
        ))

        # Remove outliers
        tmp1 <- subset(tmp, acqs == acq[1 + jj * 2 - 2])
        outl1 <- MMCD(x = tmp1$ses_01, y = tmp1$ses_02)
        tmp1 <- tmp1[outl1, ]
        tmp2 <- subset(tmp, acqs == acq[2 + jj * 2 - 2])
        outl2 <- MMCD(x = tmp2$ses_01, y = tmp2$ses_02)
        tmp2 <- tmp2[outl2, ]
        tmp <- rbind(tmp1, tmp2)

        # Here is where we do the stats for the Bland-Altman plots
        tmp$means <- rowMeans(tmp[, c(3:4)], na.rm = TRUE)
        tmp$diffs <- tmp$ses_01 - tmp$ses_02
        tmp$prct.diffs <- tmp$diffs / tmp$means * 100

        # Create empty matrices
        mean.prct.diff <- matrix(data = NA, nrow = 2, ncol = 1)
        agree.lb <- matrix(data = NA, nrow = 2, ncol = 1)
        agree.ub <- matrix(data = NA, nrow = 2, ncol = 1)
        ci.lb <- matrix(data = NA, nrow = 2, ncol = 1)
        ci.ub <- matrix(data = NA, nrow = 2, ncol = 1)

        # Calculate the mean percentage differences between scans and the 
        # corresponding limits of agreement and 95% confidence intervals
        for (ll in 1:2) {
          tmp1 <- subset(tmp, acqs == acq[ll + jj * 2 - 2])
          mean.prct.diff[ll, ] <- mean(tmp1$prct.diffs, na.rm = TRUE)
          agree.lb[ll, ] <- mean.prct.diff[ll, ] - (1.96 * sd(tmp1$prct.diffs, na.rm = TRUE))
          agree.ub[ll, ] <- mean.prct.diff[ll, ] + (1.96 * sd(tmp1$prct.diffs, na.rm = TRUE))
          ci.lb[ll, ] <- mean.prct.diff[ll, ] - sqrt(var(tmp1$prct.diffs, na.rm = TRUE) / length(tmp1))
          ci.ub[ll, ] <- mean.prct.diff[ll, ] + sqrt(var(tmp1$prct.diffs, na.rm = TRUE) / length(tmp1))
        }

        # Tibble for mean prct. differences, agreement limits, and confidence intervals
        tmp2 <- as_tibble(
          data.frame(
            acqs = acq[1:2 + (jj - 1) * 2],
            "mean.prct.diff" = mean.prct.diff,
            "agree.lb" = agree.lb,
            "agree.ub" = agree.ub,
            "ci.lb" = ci.lb,
            "ci.ub" = ci.ub
          )
        )

        # Now create the plot
        p <- ggplot(tmp, aes(x = means, y = prct.diffs, group = acqs, color = acqs)) +
          geom_point(size = 1) +
          labs(
            x = paste0("Mean of Scans 1 & 2 (", str_replace(kk, "_", "/"), ")"),
            y = "Mean percentage difference (%)"
          ) +
          theme_classic() +
          theme(
            aspect.ratio = 1,
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5, color = "black"),
            plot.margin = unit(c(0.2, 0.5, 0.2, 0.5), "cm"),
            legend.position = "bottom",
            legend.background = element_blank(),
            legend.title = element_text(face = "bold", size = 12, color = "black"),
            legend.text = element_text(size = 12, color = "black"),
            text = element_text(family = "Arial"),
            axis.ticks.length = unit(2, "mm"),
            axis.text = element_text(color = "black", size = 12),
            axis.title = element_text(color = "black", size = 14, face = "bold"),
            axis.line = element_line(color = "black", lineend = "square", linewidth = 0.25),
            axis.ticks = element_line(color = "black", linewidth = 0.25, lineend = "square")
          ) +
          ggtitle(str_replace(kk, "_", "/")) +
          scale_color_manual(
            values = c("#1b9e77", "#d95f02"),
            name = element_blank(),
            labels = c("PRESS", "sLASER")
          ) +
          scale_x_continuous(
            limits = c(
              eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[1],
              eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[2]
            ),
            breaks = seq(eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[1],
              eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[2],
              by = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[3]
            ),
            expand = c(0, 0)
          ) +
          scale_y_continuous(
            limits = c(
              eval(parse(text = paste0("yLims[[kk]]$", scheme[jj])))[1],
              eval(parse(text = paste0("yLims[[kk]]$", scheme[jj])))[2]
            ),
            breaks = seq(eval(parse(text = paste0("yLims[[kk]]$", scheme[jj])))[1],
              eval(parse(text = paste0("yLims[[kk]]$", scheme[jj])))[2],
              by = eval(parse(text = paste0("yLims[[kk]]$", scheme[jj])))[3]
            ),
            expand = c(0, 0)
          )

        p <- p + geom_segment(tmp2,
          linewidth = 0.25, lineend = "square",
          mapping = aes(
            x = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[1],
            xend = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[2],
            y = mean.prct.diff, yend = mean.prct.diff, color = acqs
          )
        ) +
          geom_segment(tmp2,
            linewidth = 0.25, linetype = 2, lineend = "square",
            mapping = aes(
              x = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[1],
              xend = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[2],
              y = agree.lb, yend = agree.lb, color = acqs
            )
          ) +
          geom_segment(tmp2,
            linewidth = 0.25, linetype = 2, lineend = "square",
            mapping = aes(
              x = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[1],
              xend = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[2],
              y = agree.ub, yend = agree.ub, color = acqs
            )
          )

        p <- p + annotate("text",
          x = eval(parse(text = paste0("xLims[[kk]]$", scheme[jj])))[2], y = mean.prct.diff,
          label = paste0(as.character(round(mean.prct.diff, 1)), "%"), size = 4, vjust = -0.2, hjust = 1,
          color = c("#1b9e77", "#d95f02"), fontface = "bold"
        )

        print(p)

        if (excl) {
          ggsave(
            file.path(out_dir, "BA_plots", paste0("BA_plot_", kk, "_", scheme[jj], "_", ii, "_excl.pdf")),
            device = cairo_pdf, height = 5, width = 5, units = "in"
          )
        } else {
          ggsave(
            file.path(out_dir, "BA_plots", paste0("BA_plot_", kk, "_", scheme[jj], "_", ii, ".pdf")),
            device = cairo_pdf, height = 5, width = 5, units = "in"
          )
        }
      }
    }
  }
}


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
