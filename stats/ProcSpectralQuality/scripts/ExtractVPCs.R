# Variance partition coefficients from lmer models ----------------------------

ExtractVPCs <- function(models, digits = 1, verbose = TRUE) {
  #' Extract VPCs from a named list of lmer models
  #'
  #' @param models  Named list of lmerMod objects (e.g., LMEM_MODELS from RunLMEM)
  #' @param digits  Number of decimal places for rounding percentages (default: 1)
  #' @param verbose If TRUE, print a formatted table per model (default: TRUE)
  #' @return A named list of data frames, one per model
  
  ExtractVPC <- function(model, digits = 1, verbose = TRUE) {
    #' Extract variance partition coefficients (VPCs) as percentages from an lmer model
    #'
    #' @param model   An lmerMod object (output of lmer() or RunLMEM())
    #' @param digits  Number of decimal places for rounding percentages (default: 1)
    #' @param verbose If TRUE, print a formatted table to the console (default: TRUE)
    #' @return A data frame with columns: grp, variance, vpc_pct
    
    vc_df     <- as.data.frame(VarCorr(model))
    groups    <- vc_df[["grp"]]
    variances <- vc_df[["vcov"]]
    vpc_pct   <- round(variances / sum(variances) * 100, digits)
    
    result <- data.frame(
      grp      = groups,
      variance = variances,
      vpc_pct  = vpc_pct,
      row.names = NULL
    )
    
    if (verbose) {
      cat("\n── Variance Partition Coefficients ──\n")
      print(result, row.names = FALSE)
      cat(sprintf("   Total variance: %.6f\n", sum(variances)))
    }
    
    invisible(result)
  }
  
  result <- lapply(names(models), function(nm) {
    if (verbose) cat(sprintf("\n── Model: %s ──", nm))
    ExtractVPC(models[[nm]], digits = digits, verbose = verbose)
  })
  
  names(result) <- names(models)
  
  return(result)
  
}
