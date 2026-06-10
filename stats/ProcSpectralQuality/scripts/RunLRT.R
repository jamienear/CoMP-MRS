# Run inference tests

RunLRT <- function(models,
                   model_contrasts,
                   run_LRTs = TRUE,
                   run_pbkrtest = FALSE,
                   verbose = FALSE) {

  if (verbose) cat("\n── Running likelihood ratio test (LRT) ──\n")

  results <- list()

  if (run_LRTs) {
    results$LRT <- mapply(
      function(large, small) {
        anova(
          models[[small]],
          models[[large]],
          test = "LRT"
        )
      },
      model_contrasts$large_models,
      model_contrasts$small_models,
      SIMPLIFY = FALSE
    )
    names(results$LRT) <- paste0(
      model_contrasts$large_models,
      "_vs_",
      model_contrasts$small_models
    )
    if (verbose) print(results$LRT)
  }

  if (run_pbkrtest) {
    results$PB <- mapply(
      function(large, small) {
        pbkrtest::PBmodcomp(
          largeModel = models[[large]],
          smallModel = models[[small]]
        )
      },
      model_contrasts$large_models,
      model_contrasts$small_models,
      SIMPLIFY = FALSE
    )
    names(results$PB) <- paste0(
      model_contrasts$small_models,
      "_vs_",
      model_contrasts$large_models
    )
    if (verbose) print(results$PB)
  }

  invisible(results)

}
