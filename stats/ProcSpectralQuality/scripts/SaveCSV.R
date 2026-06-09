### Save descriptive stats tables

SaveCSV <- function(data, out_dir) {
  
  write_csv(data$DP$LW, file.path(out_dir, "stats_DP_LW.csv"))
  write_csv(data$DP$SNR, file.path(out_dir, "stats_DP_SNR.csv"))
  write_csv(data$DP$Ratio, file.path(out_dir, "stats_DP_Ratio.csv"))
  write_csv(data$DP$Product, file.path(out_dir, "stats_DP_Product.csv"))
  
  write_csv(data$SiteID$LW, file.path(out_dir, "stats_SiteID_LW.csv"))
  write_csv(data$SiteID$SNR, file.path(out_dir, "stats_SiteID_SNR.csv"))
  write_csv(data$SiteID$Ratio, file.path(out_dir, "stats_SiteID_Ratio.csv"))
  write_csv(data$SiteID$Product, file.path(out_dir, "stats_SiteID_Product.csv"))
  
  write_csv(data$Species$LW, file.path(out_dir, "stats_Species_LW.csv"))
  write_csv(data$Species$SNR, file.path(out_dir, "stats_Species_SNR.csv"))
  write_csv(data$Species$Ratio, file.path(out_dir, "stats_Species_Ratio.csv"))
  write_csv(data$Species$Product, file.path(out_dir, "stats_Species_Product.csv"))
  
  write_csv(data$Vendor$LW, file.path(out_dir, "stats_Vendor_LW.csv"))
  write_csv(data$Vendor$SNR, file.path(out_dir, "stats_Vendor_SNR.csv"))
  write_csv(data$Vendor$Ratio, file.path(out_dir, "stats_Vendor_Ratio.csv"))
  write_csv(data$Vendor$Product, file.path(out_dir, "stats_Vendor_Product.csv"))
  
  write_csv(data$FieldStrength$LW, file.path(out_dir, "stats_FieldStrength_LW.csv"))
  write_csv(data$FieldStrength$SNR, file.path(out_dir, "stats_FieldStrength_SNR.csv"))
  write_csv(data$FieldStrength$Ratio, file.path(out_dir, "stats_FieldStrength_Ratio.csv"))
  write_csv(data$FieldStrength$Product, file.path(out_dir, "stats_FieldStrength_Product.csv"))
  
  write_csv(data$Sequence$LW, file.path(out_dir, "stats_Sequence_LW.csv"))
  write_csv(data$Sequence$SNR, file.path(out_dir, "stats_Sequence_SNR.csv"))
  write_csv(data$Sequence$Ratio, file.path(out_dir, "stats_Sequence_Ratio.csv"))
  write_csv(data$Sequence$Product, file.path(out_dir, "stats_Sequence_Product.csv"))
  
  write_csv(data$VOI$LW, file.path(out_dir, "stats_VOI_LW.csv"))
  write_csv(data$VOI$SNR, file.path(out_dir, "stats_VOI_SNR.csv"))
  write_csv(data$VOI$Ratio, file.path(out_dir, "stats_VOI_Ratio.csv"))
  write_csv(data$VOI$Product, file.path(out_dir, "stats_VOI_Product.csv"))
  
  write_csv(
    data$Species_Vendor,
    file.path(out_dir, "stats_Species_Vendor.csv")
  )
  
  write_csv(
    data$Species_VOI,
    file.path(out_dir, "stats_Species_VOI.csv")
  )
  
  write_csv(
    tibble(
      meanLW = data$all$meanLW,
      sdLW = data$all$sdLW,
      cvLW = data$all$cvLW,
      meanSNR = data$all$meanSNR,
      sdSNR = data$all$sdSNR,
      cvSNR = data$all$cvSNR,
      meanRatio = data$all$meanRatio,
      sdRatio = data$all$sdRatio,
      cvRatio = data$all$cvRatio,
      meanProduct = data$all$meanProduct,
      sdProduct = data$all$sdProduct,
      cvProduct = data$all$cvProduct
    ),
    file.path(out_dir, "stats_all.csv")
  )
  
}
