# Linear mixed-effects modeling

RunLMEM <- function(data, dv = "", rand_ef = list(), fix_ef = "") {
  
  cat("\n── Running linear mixed-effects model (LMEM) ──\n")
  
  # Fit dynamic lmer models ---------------------------------------------------
  
  fit_lmer <- function(response, random_effects, fixed_effects, data, ...) {
    #' Dynamically build and fit a linear mixed model using lmer
    
    #'
    #' @param response      Character string: name of the response variable
    #' @param fixed_effects Character vector: names of fixed-effect predictors
    #'                      (e.g., c("group", "age", "sex"))
    #' @param random_effects Named list: each element defines a random-effect term.
    #'                       The name is the grouping factor; the value is a
    #'                       character vector of random slopes (use "1" or NULL
    #'                       for intercept-only).
    #'                       Examples:
    #'                         list(subject = "1")
    #'                         list(subject = c("1", "time"), site = "1")
    #' @param data           Data frame
    #' @param ...            Additional arguments passed to lmer()
    #' @return An lmerMod object
    
    # --- validate inputs ---
    stopifnot(
      is.character(response), length(response) == 1,
      is.list(random_effects), length(random_effects) >= 1,
      !is.null(names(random_effects)),
      is.character(fixed_effects), length(fixed_effects) >= 1,
      is.data.frame(data)
    )
    
    # --- fixed-effects part ---
    if (length(fixed_effects) == 1 && nchar(fixed_effects) == 0) {
      fixed_part <- "1"
    } else {
      fixed_part <- paste(fixed_effects, collapse = " + ")
    }
    
    # --- random-effects part ---
    random_part <- vapply(
      names(random_effects),
      function(grp) {
        slopes <- random_effects[[grp]]
        if (is.null(slopes)) slopes <- "1"
        terms <- paste(slopes, collapse = " + ")
        sprintf("(%s | %s)", terms, grp)
      },
      character(1)
    )
    
    random_part <- paste(random_part, collapse = " + ")
    
    # --- assemble and fit ---
    formula_str <- paste(response, "~", fixed_part, "+", random_part)
    f <- as.formula(formula_str)
    
    message("Fitting model: ", deparse(f))
    
    m <- lme4::lmer(f, data = data, ...)
    return(m)
  }
  
  
  # Selection of optimizers to use --------------------------------------------
  
  optWrap.nloptwrap     <- lmerControl(optimizer = "nloptwrap")
  optWrap.optimx.nlminb <- lmerControl(optimizer = "optimx", optCtrl = list(method = "nlminb", eval.max = 1e5))
  optWrap.minqa.bobyqa  <- lmerControl(optimizer = "bobyqa")
  
  optimToUse <- optWrap.minqa.bobyqa
  
  
  # Remove rows with missing values -------------------------------------------
  
  data <- drop_na(data, dplyr::all_of(dv))
  
  
  # Models --------------------------------------------------------------------
  
  M <- fit_lmer(
    response = dv,
    random_effects = rand_ef,
    fixed_effects = fix_ef,
    data = data,
    REML = TRUE,
    control = optimToUse
  )
  return(M)
  
}
