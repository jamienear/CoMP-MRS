# Set in-line functions here to avoid cluttering the main analysis script

cv <- function(x, na.rm = FALSE) {
  sd(x, na.rm = na.rm) / mean(x, na.rm = na.rm)
}


# Common theme for plots ------------------------------------------------------

theme_comp <- function() {
  theme_classic() +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(
        face = "bold",
        size = 13,
        color = "black",
        family = "sans"
      ),
      axis.text = element_text(
        size = 11,
        color = "black",
        family = "sans"
      ),
      strip.text = element_text(
        face = "bold",
        size = 12,
        color = "black",
        family = "sans"
      ),
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      panel.spacing = unit(1, "lines")
    )
}


# Functions for lmer ----------------------------------------------------------

nloptFun <- function(fn, par, lower, upper, control = list(), ...) {
  for (n in names(defaultControl))
    if (is.null(control[[n]])) control[[n]] <- defaultControl[[n]]
  res <- nloptr(x0 = par, eval_f = fn, lb = lower, ub = upper, opts = control, ...)
  with(res, list(par = solution,
                 fval = objective,
                 feval = iterations,
                 conv = if (status>0) 0 else status,
                 message = message))
}

# For standardizing (z-transforming) outcome and predictor variables; aids with
# model convergence and interpretability of parameter estimates
cs. <- function(x) scale(x, center = TRUE, scale = TRUE)
