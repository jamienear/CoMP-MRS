# Export LRT results table to Word (.docx)

ExportLRTTable <- function(lrt_results,
                           out_dir  = getwd(),
                           filename = "LRT_table.docx",
                           title    = "Likelihood Ratio Tests: Linear Mixed-Effects Models") {
  #'
  #' Build a summary table from a named list of anova (LRT) objects and
  #' export it as a Word-compatible .docx file.
  #'
  #' @param lrt_results Named list of anova objects (output of RunLRT()$LRT)
  #' @param out_dir     Directory in which to save the .docx file
  #' @param filename    Output filename
  #' @param title       Table title shown above the table
  #'
  #' @return Invisibly returns the flextable object

  if (is.null(lrt_results) || length(lrt_results) == 0) {
    stop("No LRT results provided.")
  }

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # ── Build tidy data frame from list of anova objects ──────────────────────
  rows <- lapply(names(lrt_results), function(nm) {
    av  <- lrt_results[[nm]]

    # Parse model names from "large_vs_small" key
    parts       <- strsplit(nm, "_vs_")[[1]]
    large_model <- parts[1]
    small_model <- parts[2]

    # anova() returns 2 rows: row 1 = small model, row 2 = large model
    small_row <- av[1, ]
    large_row <- av[2, ]

    # Column names differ across lme4 versions:
    #   per-model df  : "npar" (>= 1.1-21) or "Df" (older)
    #   chi-square df : "Df"   (>= 1.1-21) or "Chi Df" (older)
    nparm_col <- if ("npar" %in% names(av)) "npar" else "Df"
    chidf_col <- if ("npar" %in% names(av)) "Df"   else "Chi Df"

    data.frame(
      Comparison       = nm,
      "Small model"    = small_model,
      "Large model"    = large_model,
      "df (small)"     = small_row[[nparm_col]],
      "df (large)"     = large_row[[nparm_col]],
      # "Δdf"            = large_row[[chidf_col]],
      # "AIC (small)"    = round(small_row[["AIC"]],  2),
      # "AIC (large)"    = round(large_row[["AIC"]],  2),
      # "ΔAIC"           = round(large_row[["AIC"]] - small_row[["AIC"]], 2),
      # "BIC (small)"    = round(small_row[["BIC"]],  2),
      # "BIC (large)"    = round(large_row[["BIC"]],  2),
      # "ΔBIC"           = round(large_row[["BIC"]] - small_row[["BIC"]], 2),
      "logLik (small)" = round(small_row[["logLik"]], 2),
      "logLik (large)" = round(large_row[["logLik"]], 2),
      "χ²"             = round(large_row[["Chisq"]], 3L),
      "p"              = format.pval(large_row[["Pr(>Chisq)"]], digits = 3L, eps = 0.001),
      check.names     = FALSE,
      stringsAsFactors = FALSE
    )
  })

  tbl_df <- do.call(rbind, rows)

  # ── Build flextable ────────────────────────────────────────────────────────
  ft <- flextable::flextable(tbl_df) |>
    flextable::theme_booktabs() |>
    flextable::set_caption(caption = title) |>
    flextable::add_footer_lines(
      # values = c(
      #   "ΔAIC / ΔBIC: difference between large and small model (negative = better fit for large model).",
      #   "χ² and p-value derived from likelihood ratio tests comparing nested models.",
      #   "Models fitted with REML = FALSE for valid LRT comparison."
      # )
      values = c(
        "χ² and p-value derived from likelihood ratio tests comparing nested models.",
        "Models fitted with REML = FALSE for valid LRT comparison."
      )
    ) |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::set_table_properties(layout = "autofit")

  # ── Save to Word ───────────────────────────────────────────────────────────
  out_path        <- file.path(out_dir, filename)
  sect_properties <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape")
  )
  flextable::save_as_docx(ft, path = out_path, pr_section = sect_properties)

  message("LRT table exported to: ", out_path)
  invisible(ft)
}
