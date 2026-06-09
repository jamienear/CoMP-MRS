# ============================================================
#  Outlier Exclusion in Nested Data via Median Absolute
#  Deviation (MAD) — Prepared for lmer Analysis
#
#  Strategies covered
#  ─────────────────────────────────────────────────────────
#  1. Global MAD        — single threshold across all data
#  2. Group-level MAD   — threshold computed within each
#                         level of a grouping factor
#  3. Multilevel MAD    — applied sequentially at each
#                         level of a nested hierarchy
# ============================================================
#
# Created using Claude Opus 4.6 (2026-03-30)
# Modified by Mark Mikkelsen, Ph.D. (2026-06-09)

mad_outlier_removal <- function(data,
                                outcome,
                                groups,
                                strategy   = c("multilevel", "global", "group"),
                                threshold  = 2.5,
                                thresholds = NULL, # per-level for multilevel: length(groups) + 1
                                verbose    = TRUE) {

  strategy <- match.arg(strategy)

  # ============================================================
  #  Input validation
  # ============================================================

  if (!outcome %in% names(data)) {
    stop(sprintf("Column '%s' not found in `data`.", outcome))
  }
  if (!is.character(groups) || length(groups) == 0L) {
    stop("`groups` must be a non-empty character vector of column names.")
  }
  missing_groups <- setdiff(groups, names(data))
  if (length(missing_groups) > 0L) {
    stop(sprintf("Grouping column(s) not found in `data`: %s",
                 paste(missing_groups, collapse = ", ")))
  }
  n_na <- sum(is.na(data[[outcome]]))
  if (n_na > 0L)
    warning(sprintf(
      "%d NA(s) in '%s': these rows cannot be scored and will be retained in data_clean.",
      n_na, outcome
    ))


  # ============================================================
  #  MAD core utilities
  # ============================================================

  # Compute median and scale-adjusted MAD, with mean-absolute-deviation
  # fallback when MAD is zero (e.g. >50% identical values).
  .mad_parts <- function(x) {
    med     <- stats::median(x, na.rm = TRUE)
    mad_val <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(mad_val) || mad_val == 0)
      mad_val <- mean(abs(x - med), na.rm = TRUE)
    list(med = med, mad_val = mad_val)
  }

  #' Compute the MAD-based modified Z-score for a numeric vector
  #'
  #' Formula (Iglewicz & Hoaglin 1993):
  #'   M_i = 0.6745 * (x_i - median(x)) / MAD(x)
  #'
  #' The constant 0.6745 (:= 1 / 1.4826) makes the score consistent with a
  #' normal distribution: |M_i| > threshold (default 2.5)
  #' flags an outlier.
  #'
  #' @param x         Numeric vector
  #' @param threshold Absolute M-score cut-off (default 2.5)
  #' @param constant  Consistency constant (default 0.6745)
  #' @return          Logical vector, TRUE = outlier
  mad_outlier <- function(x, threshold = 2.5, constant = 0.6745) {

    if (!is.numeric(x)) stop("`x` must be numeric.")
    if (length(x) < 3L) return(rep(FALSE, length(x)))   # too few to judge

    parts <- .mad_parts(x)
    if (!is.finite(parts$mad_val) || parts$mad_val == 0)
      return(rep(FALSE, length(x)))

    constant * abs(x - parts$med) / parts$mad_val > threshold
  }

  # Return the modified Z-score (numeric) rather than a flag
  mad_score <- function(x, constant = 0.6745) {

    if (length(x) < 3L) return(rep(NA_real_, length(x)))

    parts <- .mad_parts(x)
    if (!is.finite(parts$mad_val) || parts$mad_val == 0)
      return(rep(0, length(x)))

    constant * (x - parts$med) / parts$mad_val
  }

  # ============================================================
  #  Strategy 1: Global MAD
  # ============================================================

  #' Flag outliers using a single MAD threshold across all rows
  #'
  #' Best when the outcome is expected to be identically
  #' distributed across groups (rarely the case in nested data).
  flag_global_mad <- function(data, outcome, threshold = 2.5) {

    y   <- data[[outcome]]
    flg <- mad_outlier(y, threshold = threshold)

    data %>%
      dplyr::mutate(
        mad_score_global = mad_score(y),
        outlier_global   = flg
      )
  }

  # ============================================================
  #  Strategy 2: Group-level MAD
  # ============================================================
  
  #' Flag outliers within each level of one grouping factor
  #'
  #' Accounts for group means differing; removes observations
  #' that are extreme *relative to their own group*.
  #'
  #' @param data      Data frame
  #' @param outcome   Column name (string) of the response
  #' @param group     Column name (string) of the grouping variable
  #' @param threshold Modified Z-score cut-off
  flag_group_mad <- function(data, outcome, group, threshold = 2.5) {

    score_col  <- dplyr::sym(outcome)
    group_col  <- dplyr::sym(group)
    flag_col   <- paste0("outlier_", group)
    mscore_col <- paste0("mad_score_", group)

    data %>%
      dplyr::group_by({{ group_col }}) %>%
      dplyr::mutate(
        !!mscore_col := mad_score(!!score_col),
        !!flag_col   := mad_outlier(!!score_col, threshold = threshold)
      ) %>%
      dplyr::ungroup()
  }

  # ============================================================
  #  Strategy 3: Multilevel MAD (sequential)
  # ============================================================

  #' Apply MAD sequentially at each level of the hierarchy
  #'
  #' Step 1 — flag global extremes (rare, severe outliers).
  #' Step 2 — within each highest-level cluster, flag extremes.
  #' Step 3 — within each lowest-level cluster, flag extremes.
  #' An observation is flagged if it is extreme at *any* level.
  #'
  #' Thresholds can differ per level:
  #' E.g.:
  #'   level 1 (global)  → tighter (default 3.5)
  #'   level 2 (school)  → moderate (default 3.5)
  #'   level 3 (class)   → slightly looser (default 3.0)
  flag_multilevel_mad <- function(data,
                                  outcome,
                                  groups,     # character vector, outer → inner
                                  thresholds = NULL, # one per group + global
                                  verbose    = TRUE) {

    # Default thresholds: same for all levels
    n_levels <- length(groups) + 1   # groups + global
    if (is.null(thresholds)) thresholds <- rep(2.5, n_levels)
    stopifnot(length(thresholds) == n_levels)

    score_col <- rlang::sym(outcome)

    # Level 0: global
    data <- data %>%
      dplyr::mutate(
        mad_score_global = mad_score(!!score_col),
        .global_flag     = mad_outlier(!!score_col, threshold = thresholds[1])
      )

    # Capture per-level counts before combining (for accurate verbose output)
    level_counts <- integer(length(groups) + 1L)
    level_counts[1L] <- sum(data$.global_flag, na.rm = TRUE)

    # Level k: per group
    for (k in seq_along(groups)) {
      grp     <- rlang::sym(groups[k])
      tmp_col <- paste0(".grp_flag_", k)
      thr     <- thresholds[k + 1]

      data <- data %>%
        dplyr::group_by({{ grp }}) %>%
        dplyr::mutate(!!tmp_col := mad_outlier(!!score_col, threshold = thr)) %>%
        dplyr::ungroup()

      level_counts[k + 1L] <- sum(data[[tmp_col]], na.rm = TRUE)
    }

    # Combine: flagged at ANY level
    flag_cols <- c(".global_flag",
                   paste0(".grp_flag_", seq_along(groups)))

    data <- data %>%
      dplyr::mutate(
        outlier_multilevel = rowSums(dplyr::across(dplyr::all_of(flag_cols)),
                                     na.rm = TRUE) > 0
      ) %>%
      dplyr::select(-dplyr::all_of(flag_cols))

    if (verbose) {
      n_flagged <- sum(data$outlier_multilevel)
      pct       <- round(100 * n_flagged / nrow(data), 2)
      cat(sprintf("Multilevel MAD: %d / %d flagged (%.1f%%)\n",
                  n_flagged, nrow(data), pct))

      # Per-level breakdown (captured before temp columns were dropped)
      cat("  Global level:", level_counts[1L], "\n")
      for (k in seq_along(groups)) {
        cat(sprintf("  %s level: %d\n", groups[k], level_counts[k + 1L]))
      }
      cat("\n")
    }

    data
  }

  # ============================================================
  # Run outlier removal using chosen strategy
  # ============================================================

  if (verbose) {
    cat(sprintf("\n─── outlier_removal | strategy = '%s' ───\n",
                strategy))
    cat(sprintf("   Input rows: %d | outcome: '%s'\n", nrow(data), outcome))
  }

  # Apply chosen strategy
  if (strategy == "global") {
    data <- flag_global_mad(data, outcome, threshold)
    flag_col <- "outlier_global"

  } else if (strategy == "group") {
    if (length(groups) > 1L)
      warning(sprintf(
        "`strategy = 'group'` only uses the first group ('%s'); remaining groups ignored.",
        groups[1]
      ))
    data     <- flag_group_mad(data, outcome, groups[1], threshold)
    flag_col <- paste0("outlier_", groups[1])

  } else if (strategy == "multilevel") {
    thresholds <- if (is.null(thresholds))
      rep(threshold, length(groups) + 1) else thresholds
    if (length(thresholds) != length(groups) + 1L)
      stop(sprintf(
        "`thresholds` must have length %d (global + one per group).",
        length(groups) + 1L
      ))
    data <- flag_multilevel_mad(data, outcome, groups,
                                thresholds = thresholds,
                                verbose    = verbose)
    flag_col <- "outlier_multilevel"
  }

  # NA flags (from NA outcomes) are treated as not-flagged: row is retained
  flag_vec   <- data[[flag_col]]
  data_clean <- dplyr::filter(data, !flag_vec | is.na(flag_vec))

  # Build report
  report <- list(
    strategy     = strategy,
    threshold    = threshold,
    thresholds   = if (strategy == "multilevel") thresholds else NULL,
    n_original   = nrow(data),
    n_flagged    = sum(data[[flag_col]], na.rm = TRUE),
    n_clean      = nrow(data_clean),
    pct_removed  = round(100 * mean(data[[flag_col]], na.rm = TRUE), 2),
    flag_column  = flag_col
  )

  if (verbose) {
    cat(sprintf("   Flagged: %d (%.1f%%) | Retained: %d\n\n",
                report$n_flagged, report$pct_removed, report$n_clean))
  }

  list(
    data         = data,
    data_clean   = data_clean,
    outlier_flag = data[[flag_col]],
    report       = report
  )

}
