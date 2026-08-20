# Export a modelsummary table of LMEM results to Word (.docx)

# ── Custom glance method ─────────────────────────────────────────────────────
# modelsummary dispatches glance_custom(model) and merges the result into the
# GOF section. Defining this for lmerMod injects conditional R² (fixed +
# random effects) so it can be referenced in gof_map as "r2.conditional".
# lm objects (e.g., null model) fall back to the default empty data frame,
# leaving the R² cell blank for that column.
glance_custom.lmerMod <- function(x, ...) {
  r2 <- performance::r2(x)
  data.frame(r2.conditional = as.numeric(r2$R2_conditional))
}

ExportModelTable_old <- function(models,
                                 vpcs = NULL,
                                 out_dir,
                                 filename = "LMEM_table_old.docx",
                                 title = "Linear Mixed-Effects Models: Norm. SNR/LW") {
  #'
  #' Build a modelsummary comparison table from a named list of lm/lmerMod
  #' objects and export it as a Word-compatible .docx file.
  #'
  #' @param models   Named list of model objects (lm, lmerMod, etc.)
  #' @param vpcs     Named list of VPC data frames from ExtractVPCs() (optional).
  #'                 Names must match names(models). Models absent from vpcs
  #'                 (e.g., lm null model) receive empty cells in VPC rows.
  #' @param out_dir  Directory in which to save the .docx file
  #' @param filename Output filename (default: "LMEM_table.docx")
  #' @param title    Table title shown above the table
  #'
  #' @return Invisibly returns the flextable object
  
  if (is.null(models) || length(models) == 0) {
    stop("No models provided. Please supply a named list of model objects.")
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # ── Goodness-of-fit rows to include ────────────────────────────────────────
  gof_map <- dplyr::tribble(
    ~raw,             ~clean,                  ~fmt,
    # "nobs",           "Observations",          "%.0f",
    "r2.conditional", "R\u00B2 (conditional)", "%.3f"
  )
  
  # ── Build VPC add_rows ─────────────────────────────────────────────────────
  # add_rows must have 1 label column + 1 column per model (positional match).
  # Empty string ("") is used where a grouping factor does not appear in a model.
  vpc_rows <- NULL
  if (!is.null(vpcs) && length(vpcs) > 0) {
    n_models <- length(models)
    model_names <- names(models)
    
    # All unique grouping factors across all VPC-containing models, preserving
    # the order they first appear (Residual will naturally come last)
    all_grps <- unique(unlist(lapply(vpcs, function(df) df$grp)))
    
    rows <- lapply(all_grps, function(grp) {
      vals <- vapply(
        model_names,
        function(nm) {
          vpc_df <- vpcs[[nm]]
          if (is.null(vpc_df)) return("")
          idx <- match(grp, vpc_df$grp)
          if (is.na(idx)) return("")
          sprintf("%.1f%%", vpc_df$vpc_pct[idx])
        },
        character(1)
      )
      c(paste0("VPC: ", grp), unname(vals))
    })
    
    vpc_rows <- as.data.frame(
      do.call(rbind, rows),
      stringsAsFactors = FALSE,
      check.names      = FALSE
    )
    # Use positional column names ("1", "2", ...) to match modelsummary's
    # internal column indexing rather than relying on display header names
    colnames(vpc_rows) <- c("term", as.character(seq_len(n_models)))
  }
  
  # ── Notes ──────────────────────────────────────────────────────────────────
  notes <- c(
    "Standard errors in parentheses.",
    "R\u00B2 (conditional) accounts for both fixed and random effects (performance::r2)."
  )
  if (!is.null(vpc_rows)) {
    notes <- c(
      notes,
      "VPC = Variance partition coefficient (% of total variance attributed to each grouping factor)."
    )
  }
  
  # ── Build table ────────────────────────────────────────────────────────────
  tbl <- modelsummary(
    models,
    output    = "flextable",
    # statistic = "({std.error})",
    statistic = "p.value",
    gof_map   = gof_map,
    add_rows  = vpc_rows,
    align     = str_c(c("l", rep("l", length(models))), collapse = ""),
    title     = title,
    notes     = notes
  )
  
  # ── Light formatting ───────────────────────────────────────────────────────
  tbl <- tbl |>
    flextable::fontsize(size = 10, part = "all") |>
    flextable::font(fontname = "Arial", part = "all") |>
    flextable::set_table_properties(layout = "autofit")
  
  # ── Save to Word ───────────────────────────────────────────────────────────
  out_path <- file.path(out_dir, filename)
  sect_properties <- officer::prop_section(
    page_size = officer::page_size(orient = "landscape")
  )
  flextable::save_as_docx(tbl, path = out_path, pr_section = sect_properties)
  
  message("Model table exported to: ", out_path)
  invisible(tbl)
}
