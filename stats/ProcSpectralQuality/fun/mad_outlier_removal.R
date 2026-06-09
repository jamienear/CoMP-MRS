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
                                strategy  = "multilevel",
                                threshold = 2.5,
                                verbose   = TRUE) {
  
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
  
  # ============================================================
  #  MAD core utilities
  # ============================================================
  
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
    
    if (!is.numeric(x)) {
      stop("`x` must be numeric.")
    }
    
    if (length(x) < 3L) return(rep(FALSE, length(x)))   # too few to judge
    
    med     <- stats::median(x, na.rm = TRUE)
    mad_val <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
    
    # Fallback to mean absolute deviation if raw MAD is zero/non-finite
    if (!is.finite(mad_val) || mad_val == 0) {
      mad_val <- mean(abs(x - med), na.rm = TRUE)
    }
    if (!is.finite(mad_val) || mad_val == 0) return(rep(FALSE, length(x)))
    
    m_score <- constant * abs(x - med) / mad_val
    outlier <- m_score > threshold
    outlier
  }
  
  # Return the modified Z-score (numeric) rather than a flag
  mad_score <- function(x, constant = 0.6745) {
    
    if (length(x) < 3L) return(rep(NA_real_, length(x)))
    
    med     <- stats::median(x, na.rm = TRUE)
    mad_val <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)
    if (!is.finite(mad_val) || mad_val == 0)
      mad_val <- mean(abs(x - med), na.rm = TRUE)
    if (!is.finite(mad_val) || mad_val == 0) return(rep(0, length(x)))
    
    score <- constant * (x - med) / mad_val
    score
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
    
    score_sym <- rlang::sym(outcome)
    
    # Level 0: global
    data <- data %>%
      dplyr::mutate(
        .global_flag = mad_outlier(!!score_sym, threshold = thresholds[1])
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
        dplyr::mutate(!!tmp_col := mad_outlier(!!score_sym, threshold = thr)) %>%
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
  
  stopifnot(strategy %in% c("global", "group", "multilevel"))
  
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
    # Apply at the outermost grouping level only
    data     <- flag_group_mad(data, outcome, groups[1], threshold)
    flag_col <- paste0("outlier_", groups[1])
    
  } else if (strategy == "multilevel") {
    data <- flag_multilevel_mad(data, outcome, groups,
                                thresholds = rep(threshold,
                                                 length(groups) + 1),
                                verbose    = verbose)
    flag_col <- "outlier_multilevel"
  }
  
  data_clean <- dplyr::filter(data, !.data[[flag_col]])
  
  # Build report
  report <- list(
    strategy     = strategy,
    threshold    = threshold,
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
