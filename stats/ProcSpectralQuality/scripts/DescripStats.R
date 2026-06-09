# Descriptive statistics for each grouping variable and their combinations.

DescripStats <- function(data) {
  
  STATS <- list()
  
  ### Descriptive stats function ----------------------------------------------
  
  make_stats <- function(data, group_var, value_var, prefix) {
    data %>%
      dplyr::group_by(across(all_of(group_var))) %>%
      dplyr::summarise(
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
  
  ### Species -----------------------------------------------------------
  
  STATS$Species <- list(
    LW      = make_stats(data, "Species", "LW_norm", "LW"),
    SNR     = make_stats(data, "Species", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Species", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Species", "SNR_LW_Product_norm", "Product")
  )
  
  ### Vendor ----------------------------------------------------------------
  
  STATS$Vendor <- list(
    LW      = make_stats(data, "Vendor", "LW_norm", "LW"),
    SNR     = make_stats(data, "Vendor", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Vendor", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Vendor", "SNR_LW_Product_norm", "Product")
  )
  
  ### FieldStrength -----------------------------------------------------------------
  
  STATS$FieldStrength <- list(
    LW      = make_stats(data, "FieldStrength", "LW_norm", "LW"),
    SNR     = make_stats(data, "FieldStrength", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "FieldStrength", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "FieldStrength", "SNR_LW_Product_norm", "Product")
  )
  
  ### Sequence --------------------------------------------------------------
  
  STATS$Sequence <- list(
    LW      = make_stats(data, "Sequence", "LW_norm", "LW"),
    SNR     = make_stats(data, "Sequence", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "Sequence", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "Sequence", "SNR_LW_Product_norm", "Product")
  )
  
  ### VOI -----------------------------------------------------------
  
  STATS$VOI <- list(
    LW      = make_stats(data, "VOI", "LW_norm", "LW"),
    SNR     = make_stats(data, "VOI", "SNR_norm", "SNR"),
    Ratio   = make_stats(data, "VOI", "SNR_LW_Ratio_norm", "Ratio"),
    Product = make_stats(data, "VOI", "SNR_LW_Product_norm", "Product")
  )
  
  ### Species x Vendor ------------------------------------------------
  
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
      .groups = "drop"
    )
  
  ### Species x VOI -------------------------------------------
  
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
      .groups = "drop"
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
