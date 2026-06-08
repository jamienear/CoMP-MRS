# Export a modelsummary table of LMEM results to Word (.docx)

# ── Custom glance method ─────────────────────────────────────────────────────
# modelsummary dispatches glance_custom(model) and merges the result into the
# GOF section. Defining this for lmerMod injects conditional R² (fixed +
# random effects) so it can be referenced in gof_map as "r2.conditional".
# lm objects (e.g., null model) fall back to the default empty data frame,
# leaving the R² cell blank for that column.

# Export a modelsummary table of LMEM results to Word (.docx)

ExportModelTable <- function(models,
                             vpcs = NULL,
                             out_dir = getwd(),
                             filename = "LMEM_table.docx",
                             title = flextable::as_paragraph(
                               "Linear Mixed-Effects Models: SNR/LW",
                               flextable::as_sub("norm"))) {
  #'
  #' Build a modelsummary comparison table from a named list of lmerMod
  #' objects and export it as a Word-compatible .docx file.
  #'
  #' @param models   Named list of model objects
  #' @param vpcs     Named list of VPC data frames from ExtractVPCs() (optional)
  #' @param out_dir  Directory in which to save the .docx file
  #' @param filename Output filename
  #' @param title    Table title shown above the table
  #'
  #' @return Invisibly returns the flextable object
  
  if (is.null(models) || length(models) == 0) {
    stop("No models provided. Please supply a named list of model objects.")
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # ── Extract conditional R² for each model ─────────────────────────────────
  r2_conditional <- vapply(
    models,
    function(m) {
      out <- tryCatch(
        performance::r2(m)$R2_conditional,
        error = function(e) NA_real_
      )
      if (length(out) == 0 || is.null(out)) {
        NA_real_
      } else {
        as.numeric(out[[1]])
      }
    },
    numeric(1)
  )
  
  # ── Determine VPC columns ──────────────────────────────────────────────────
  all_grps      <- NULL
  vpc_col_names <- NULL
  if (!is.null(vpcs) && length(vpcs) > 0) {
    all_grps      <- unique(unlist(lapply(vpcs, function(df) df$grp)))
    all_grps      <- c("Residual", setdiff(all_grps, "Residual"))
    vpc_col_names <- paste0("VPC: ", all_grps)
  }
  
  # ── Notes ──────────────────────────────────────────────────────────────────
  notes <- c(
    "Standard errors in parentheses.",
    "R\u00B2 (conditional) accounts for both fixed and random effects."
  )
  
  if (!is.null(vpc_col_names)) {
    notes <- c(
      notes,
      "VPC = Variance partition coefficient (% of total variance attributed to each grouping factor)."
    )
  }
  
  # ── Coefficient ordering ───────────────────────────────────────────────────
  all_coefs <- unique(unlist(lapply(models, function(m) modelsummary::get_estimates(m)$term)))
  priority  <- c("(Intercept)", "SD (Observations)")
  coef_map  <- setNames(
    c(priority, setdiff(all_coefs, priority)),
    c(priority, setdiff(all_coefs, priority))
  )
  
  # ── Build base table WITHOUT GOF rows ──────────────────────────────────────
  tbl <- modelsummary::modelsummary(
    models,
    output    = "flextable",
    statistic = "p.value",
    coef_map  = coef_map,
    shape     = model + term ~ statistic,
    title     = NULL,
    notes     = NULL
  )
  
  tbl <- flextable::set_caption(tbl, caption = title)
  
  # ── Inject R² and VPC columns into body ────────────────────────────────────
  bd <- tbl$body$data
  
  # First visible column is the model column
  orig_model_col <- bd[[tbl$col_keys[[1L]]]]
  model_col <- orig_model_col
  
  # Fill down blank model labels for lookup purposes
  for (i in seq_len(length(model_col))[-1L]) {
    if (nchar(trimws(model_col[[i]])) == 0L) {
      model_col[[i]] <- model_col[[i - 1L]]
    }
  }
  
  # Add conditional R² column after p
  # Show it only on the first visible row of each model block
  bd[["R² (conditional)"]] <- vapply(
    seq_along(model_col),
    function(i) {
      if (nchar(trimws(orig_model_col[[i]])) == 0L) return("")
      mdl <- model_col[[i]]
      val <- r2_conditional[[mdl]]
      if (is.na(val)) return("")
      sprintf("%.3f", val)
    },
    character(1)
  )
  
  # Add VPC columns
  if (!is.null(vpc_col_names)) {
    for (grp in all_grps) {
      col_nm <- paste0("VPC: ", grp)
      bd[[col_nm]] <- vapply(
        seq_along(model_col),
        function(i) {
          if (nchar(trimws(orig_model_col[[i]])) == 0L) return("")
          vpc_df <- vpcs[[model_col[[i]]]]
          if (is.null(vpc_df)) return("")
          idx <- match(grp, vpc_df$grp)
          if (is.na(idx)) return("")
          sprintf("%.1f%%", vpc_df$vpc_pct[idx])
        },
        character(1)
      )
    }
  }
  
  # Rebuild table with desired column order:
  # model | term | Est. | p | R² (conditional) | VPC...
  new_col_keys <- c(tbl$col_keys, "R² (conditional)")
  if (!is.null(vpc_col_names)) {
    new_col_keys <- c(new_col_keys, vpc_col_names)
  }
  
  tbl <- flextable::flextable(
    bd,
    col_keys = new_col_keys
  ) |>
    flextable::theme_booktabs() |>
    flextable::set_caption(caption = title) |>
    flextable::add_footer_lines(values = notes)
  
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