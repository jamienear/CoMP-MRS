# Descriptive statistics for each grouping variable and their combinations.

DescripStats <- function(data) {
  
  STATS <- list()
  
  ### Descriptive stats function ----------------------------------------------
  
  make_stats <- function(data, group_var, value_var, prefix) {
    data %>%
      dplyr::group_by(
        across(all_of(group_var))
      ) %>%
      dplyr::summarise(
        mean    = mean(.data[[value_var]], na.rm = TRUE),
        sd      = sd(.data[[value_var]], na.rm = TRUE),
        cv      = cv(.data[[value_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      rename(
        !!paste0("mean", prefix) := mean,
        !!paste0("sd", prefix)   := sd,
        !!paste0("cv", prefix)   := cv
      )
  }

  # The percent difference is calculated using the median because some of the  
  # variables for which we calculate SNR/LW_norm are NOT normally distributed,  
  # although we assume their normality for the LMEM
  # This function only works for 2-level variables!
  percent_difference <- function(data, group_var, value_var, lvl1, lvl2) {
    m1 <- median(data[[value_var]][data[[group_var]] == lvl1], na.rm = TRUE)
    m2 <- median(data[[value_var]][data[[group_var]] == lvl2], na.rm = TRUE)
    (m1 - m2) / m2 * 100
  }

  ### DP ----------------------------------------------------------------------
  
  STATS$DP <- list(
    LW      = make_stats(data, "DP", "LW_norm", "LW"),
    SNR     = make_stats(data, "DP", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "DP", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "DP", "SNR_LW_Product_norm", "Product")
  )
  
  ### SiteID ------------------------------------------------------------------
  
  STATS$SiteID <- list(
    LW      = make_stats(data, "SiteID", "LW_norm", "LW"),
    SNR     = make_stats(data, "SiteID", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "SiteID", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "SiteID", "SNR_LW_Product_norm", "Product")
  )
  
  ### Species -----------------------------------------------------------------
  
  STATS$Species <- list(
    LW      = make_stats(data, "Species", "LW_norm", "LW"),
    SNR     = make_stats(data, "Species", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Species", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Species", "SNR_LW_Product_norm", "Product")
  )
  
  ### Vendor ------------------------------------------------------------------
  
  STATS$Vendor <- list(
    LW      = make_stats(data, "Vendor", "LW_norm", "LW"),
    SNR     = make_stats(data, "Vendor", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Vendor", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Vendor", "SNR_LW_Product_norm", "Product")
  )
  
  ### FieldStrength -----------------------------------------------------------
  
  STATS$FieldStrength <- list(
    LW      = make_stats(data, "FieldStrength", "LW_norm", "LW"),
    SNR     = make_stats(data, "FieldStrength", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "FieldStrength", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "FieldStrength", "SNR_LW_Product_norm", "Product")
  )
  
  ### Sequence ----------------------------------------------------------------
  
  STATS$Sequence <- list(
    LW      = make_stats(data, "Sequence", "LW_norm", "LW"),
    SNR     = make_stats(data, "Sequence", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Sequence", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Sequence", "SNR_LW_Product_norm", "Product")
  )
    
  ### Sequence percentage difference (LASER vs. PRESS) -----------------------

  STATS$Sequence_percent_difference$LvP <- list(
    LW      = percent_difference(data, "Sequence", "LW_norm", "LASER", "PRESS"),
    SNR     = percent_difference(data, "Sequence", "SNR_norm", "LASER", "PRESS"),
    Ratio   = percent_difference(data, "Sequence", "SNR_LW_Ratio_norm", "LASER", "PRESS"),
    Product = percent_difference(data, "Sequence", "SNR_LW_Product_norm", "LASER", "PRESS")
  )
     
  ### Sequence percentage difference (LASER vs. STEAM) -----------------------

  STATS$Sequence_percent_difference$LvSt <- list(
    LW      = percent_difference(data, "Sequence", "LW_norm", "LASER", "STEAM"),
    SNR     = percent_difference(data, "Sequence", "SNR_norm", "LASER", "STEAM"),
    Ratio   = percent_difference(data, "Sequence", "SNR_LW_Ratio_norm", "LASER", "STEAM"),
    Product = percent_difference(data, "Sequence", "SNR_LW_Product_norm", "LASER", "STEAM")
  ) 

  ### Sequence percentage difference (LASER vs. SPECIAL) -----------------------

  STATS$Sequence_percent_difference$LvSp <- list(
    LW      = percent_difference(data, "Sequence", "LW_norm", "LASER", "SPECIAL"),
    SNR     = percent_difference(data, "Sequence", "SNR_norm", "LASER", "SPECIAL"),
    Ratio   = percent_difference(data, "Sequence", "SNR_LW_Ratio_norm", "LASER", "SPECIAL"),
    Product = percent_difference(data, "Sequence", "SNR_LW_Product_norm", "LASER", "SPECIAL")
  )
  
  ### VOI ---------------------------------------------------------------------
  
  STATS$VOI <- list(
    LW      = make_stats(data, "VOI", "LW_norm", "LW"),
    SNR     = make_stats(data, "VOI", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "VOI", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "VOI", "SNR_LW_Product_norm", "Product")
  )
  
  ### Cryoprobe ---------------------------------------------------------------
  
  STATS$Cryoprobe <- list(
    LW      = make_stats(data, "Cryoprobe", "LW_norm", "LW"),
    SNR     = make_stats(data, "Cryoprobe", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Cryoprobe", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Cryoprobe", "SNR_LW_Product_norm", "Product")
  )

  ### Cryoprobe percentage difference (TRUE vs FALSE) -----------------------------

  STATS$Cryoprobe_percent_difference <- list(
    LW      = percent_difference(data, "Cryoprobe", "LW_norm", TRUE, FALSE),
    SNR     = percent_difference(data, "Cryoprobe", "SNR_norm", TRUE, FALSE),
    Ratio   = percent_difference(data, "Cryoprobe", "SNR_LW_Ratio_norm", TRUE, FALSE),
    Product = percent_difference(data, "Cryoprobe", "SNR_LW_Product_norm", TRUE, FALSE)
  )
  
  ### ShimMethod ---------------------------------------------------------------
  
  STATS$ShimMethod <- list(
    LW      = make_stats(data, "ShimMethod", "LW_norm", "LW"),
    SNR     = make_stats(data, "ShimMethod", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "ShimMethod", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "ShimMethod", "SNR_LW_Product_norm", "Product")
  )
  
  ### ShimMethod percentage difference (FAST(EST)MAP vs MAPSHIM) -----------------------------
  
  STATS$ShimMethod_percent_difference <- list(
    LW      = percent_difference(data, "ShimMethod", "LW_norm", "MAPSHIM", "FAST(EST)MAP"),
    SNR     = percent_difference(data, "ShimMethod", "SNR_norm", "MAPSHIM", "FAST(EST)MAP"),
    Ratio   = percent_difference(data, "ShimMethod", "SNR_LW_Ratio_norm", "MAPSHIM", "FAST(EST)MAP"),
    Product = percent_difference(data, "ShimMethod", "SNR_LW_Product_norm", "MAPSHIM", "FAST(EST)MAP")
  )
  
  ### Species x Vendor --------------------------------------------------------
  
  STATS$Species_Vendor <- data %>%
    dplyr::group_by(Vendor, Species) %>%
    dplyr::summarise(
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
      .groups     = "drop"
    )
  
  ### Species x VOI -----------------------------------------------------------
  
  STATS$Species_VOI <- data %>%
    dplyr::group_by(Species, VOI) %>%
    dplyr::summarise(
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
      .groups     = "drop"
    )
  
  ### All ---------------------------------------------------------------------
  
  STATS$all <- list(
    meanLW      = mean(data$LW_norm, na.rm = TRUE),
    sdLW        = sd(data$LW_norm, na.rm = TRUE),
    cvLW        = cv(data$LW_norm, na.rm = TRUE),
    meanSNR     = mean(data$SNR_norm, na.rm = TRUE),
    sdSNR       = sd(data$SNR_norm, na.rm = TRUE),
    cvSNR       = cv(data$SNR_norm, na.rm = TRUE),
    meanRatio   = mean(data$SNR_LW_Ratio_norm, na.rm = TRUE),
    sdRatio     = sd(data$SNR_LW_Ratio_norm, na.rm = TRUE),
    cvRatio     = cv(data$SNR_LW_Ratio_norm, na.rm = TRUE),
    meanProduct = mean(data$SNR_LW_Product_norm, na.rm = TRUE),
    sdProduct   = sd(data$SNR_LW_Product_norm, na.rm = TRUE),
    cvProduct   = cv(data$SNR_LW_Product_norm, na.rm = TRUE)
  )
  
  return(STATS)
  
}
