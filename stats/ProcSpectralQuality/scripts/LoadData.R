# Load CoMP-MRS data and perform initial cleaning and outlier removal

LoadData <- function(csv_file = "CoMP_MRS_Rstats_input.csv",
                     outl_rm_strategy = "group",
                     verbose = TRUE
) {
  
  # Load the participants.csv file with appropriate column types --------------
  
  DATA <- read_csv(
    file = csv_file,
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
      LW = col_double(),
      SNR = col_double(),
      SNR_LW_Ratio = col_double(),
      CompCheck = col_factor(),
      MRSshim = col_factor(),
      Cryoprobe = col_logical()
    )
  )

  # Change variable names -----------------------------------------------------
  
  DATA <- DATA %>%
    rename(
      Vendor = MRvendor,
      FieldStrength_orig = MRfield,
      Averages = MRaverages,
      ShimMethod = MRSshim,
      Sequence_orig = MRsequence,
      VOI_orig = MRbrainregion,
      VoxelVolume = MRvoxelvolume,
      Species = AnimalSpecies,
      Sex = AnimalSex,
      Age = AnimalAge,
      Strain = AnimalStrain,
      Weight = AnimalWeight,
      SoftwareVer = MRsoftwareversion
    )
  
  # Remove DP01 (example only) ------------------------------------------------
  
  DATA <- DATA %>%
    filter(DP != "DP01") 
  
  # Rename some data entries --------------------------------------------------
  
  # Remove "compMR" prefix from CompID to make it more concise (optional)
  # One subject has "Other" as a shim method; change to "MAPSHIM" for
  # consistency with others in the DP
  # Change "FASTMAP-FASTESTMAP" to "FAST(EST)MAP" for consistency with others
  # in the DP
  DATA <- DATA %>%
    mutate(
      CompID = as.factor(str_remove_all(CompID, "compMR"))
    ) %>%
    mutate(
      ShimMethod = as.factor(if_else(ShimMethod == "Other", "MAPSHIM", ShimMethod))
    ) %>%
    mutate(
      ShimMethod = as.factor(if_else(ShimMethod == "FASTMAP-FASTESTMAP", "FAST(EST)MAP", ShimMethod))
    )
  
  # Reorder all factors alphabetically / numerically ascending ----------------
  
  DATA <- DATA %>%
    mutate(across(where(is.factor), ~ fct_relevel(., sort(levels(.))))) %>%
    mutate(
      VOI_orig = fct_relevel(
        VOI_orig,
        "Lhippocampus",
        "Rhippocampus",
        "Lstriatum",
        "Rstriatum"
      )
    ) %>%
    mutate(
      Sequence_orig = fct_relevel(
        Sequence_orig,
        "LASER",
        "sLASER",
        "PRESS",
        "SPECIAL",
        "sSPECIAL",
        "STEAM"
      )
    ) %>%
    mutate(Cryoprobe = factor(Cryoprobe, levels = c(TRUE, FALSE)))


  # Create new variables where some levels are collapsed ----------------------

  lookup_seq <- tibble::tribble(
    ~from,       ~to,
    "sLASER",   "LASER",
    "sSPECIAL", "SPECIAL"
  )

  lookup_voi <- tibble::tribble(
    ~from,            ~to,
    "Lhippocampus",   "Hippocampus",
    "Rhippocampus",   "Hippocampus",
    "Lstriatum",      "Striatum",
    "Rstriatum",      "Striatum"
  )

  DATA <- DATA %>%
    mutate(
      Sequence = factor(recode_values(
        Sequence_orig,
        from = lookup_seq$from,
        to = lookup_seq$to,
        default = Sequence_orig
      )),
      VOI = factor(recode_values(
        VOI_orig,
        from = lookup_voi$from,
        to = lookup_voi$to,
        default = VOI_orig
      )),
      # FieldStrength = factor(case_when(
      #   FieldStrength_orig %in% c(11.7, 14.1, 15.2) ~ "11.7+",
      #   TRUE ~ as.character(FieldStrength_orig)
      # )),
      FieldStrength = case_when(
        FieldStrength_orig %in% c(14.1, 15.2) ~ 11.7,
        TRUE ~ as.numeric(FieldStrength_orig)
      )
    ) %>%
    relocate(Sequence, .after = Sequence_orig) %>%
    relocate(VOI, .after = VOI_orig) %>%
    relocate(FieldStrength, .after = FieldStrength_orig)
  
  # Spectral quality metrics and normalization --------------------------------
  
  DATA <- DATA %>%
    mutate(
      SNR_LW_Ratio_norm = SNR_LW_Ratio / (sqrt(Averages) * VoxelVolume),
      SNR_LW_Product_norm = (LW * SNR) / (sqrt(Averages) * VoxelVolume),
      LW_norm = LW / (sqrt(Averages) * VoxelVolume),
      SNR_norm = SNR / (sqrt(Averages) * VoxelVolume)
    )

  # Clean up data -------------------------------------------------------------

  # Ensure that all missing values are represented as NA (in case there are any
  # other representations of missing data in the original CSV)
  DATA[is.na(DATA)] <- NA

  ### Remove datasets that failed QC ------------------------------------------

  DATA <- DATA %>%
    filter(CompCheck != "fail") %>%
    select(-CompCheck) # Remove the CompCheck column as it's no longer needed

  ### Remove outliers based on the MAD method ---------------------------------
  # This will also remove any rows where LW or SNR is NA, which is appropriate for outlier removal

  DATA_clean <- mad_outlier_removal(
    data = DATA,
    outcome = "SNR_LW_Ratio_norm",
    groups = "DP", # outer → inner
    strategy = outl_rm_strategy, # options: "global", "group", "multilevel"
    threshold = 2.5,
    verbose = verbose
  )
  DATA_orig <- DATA
  DATA <- DATA_clean$data_clean # Update DATA to the cleaned dataset

  # Optional but recommended:
  # aggregate to one row per DP for cleaner plots
  DATA_DP <- DATA %>%
    dplyr::group_by(
      DP,
      SiteID,
      Species,
      Sex,
      Vendor,
      FieldStrength_orig,
      FieldStrength,
      Sequence_orig,
      Sequence,
      VOI_orig,
      VOI
    ) %>%
    dplyr::summarise(
      LW_norm = mean(LW_norm, na.rm = TRUE),
      SNR_norm = mean(SNR_norm, na.rm = TRUE),
      SNR_LW_Product_norm = mean(SNR_LW_Product_norm, na.rm = TRUE),
      SNR_LW_Ratio_norm = mean(SNR_LW_Ratio_norm, na.rm = TRUE),
      .groups = "drop"
    )

  return(
    list(
      data = DATA,
      data_orig = DATA_orig,
      data_by_dp = DATA_DP
    )
  )
  
}
