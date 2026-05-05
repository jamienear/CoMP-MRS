# Box plots faceted by a third variable

PlotFacetBoxPlots <- function(data, out_dir, x_var, y_vars, facet_var) {

  # Normalize x_var / facet_var to lists of lists with $var (and optional $label)
  normalize_var_list <- function(v) {
    if (is.character(v)) {
      lapply(v, function(s) list(var = s, label = s))
    } else {
      lapply(v, function(item) {
        if (is.null(item$label)) item$label <- item$var
        item
      })
    }
  }

  x_vars    <- normalize_var_list(x_var)
  facet_var <- normalize_var_list(facet_var)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  } else {
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir)
  }

  # Faceted box plot function --------------------------------------------------

  make_facet_box_plot <- function(
    data, x_var, x_label, y_var, y_label,
    facet_var, facet_label, out_dir, file_name
  ) {

    every_nth <- function(n) {
      function(x) x[seq(1, length(x), by = n)]
    }

    plot_data <- data %>%
      dplyr::filter(
        !is.na(.data[[x_var]]),
        !is.na(.data[[y_var]]),
        !is.na(.data[[facet_var]])
      ) %>%
      dplyr::mutate(
        !!x_var := factor(
          .data[[x_var]], levels = sort(unique(.data[[x_var]]))
        ),
        !!facet_var := factor(
          .data[[facet_var]], levels = sort(unique(.data[[facet_var]]))
        )
      )

    n_levels <- nlevels(plot_data[[x_var]])
    nth <- max(1, round(n_levels / 50))

    if (n_levels > 10) {
      a <- 45
      j <- 1
    } else {
      a <- 0
      j <- 0.5
    }

    n_facets <- nlevels(plot_data[[facet_var]])
    plot_width <- max(8, 1.5 * n_facets)

    p <- ggplot(
      plot_data,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]],
      )
    ) +
      geom_boxplot(
        fill = "#4C77C2",
        width = 0.65,
        alpha = 0.9,
        outlier.shape = 21,
        outlier.size = 1.8,
        outlier.stroke = 0.3,
        color = "black",
        linewidth = 0.4
      ) +
      geom_point(
        alpha = 0.85,
        color = "black",
        shape = 16
      ) +
      facet_wrap(reformulate(facet_var), scales = "free_x") +
      labs(
        title = paste(y_label, "by", x_label, "and", facet_label),
        x = x_label,
        y = y_label
      ) +
      scale_x_discrete(
        breaks = every_nth(nth)
      ) +
      theme_comp() +
      theme(
        axis.text.x = element_text(angle = a, hjust = j)
      )

    print(p)

    ggsave(
      filename = file.path(out_dir, file_name),
      plot = p,
      width = plot_width,
      height = 5,
      units = "in",
      dpi = 300
    )
  }


  # Make faceted box plots -----------------------------------------------------

  all_plots <- list()

  for (f in facet_var) {
    for (x in x_vars) {
      for (y in y_vars) {
        file_name <- paste0(
          "facet_boxplot_", y$var, "_by_", x$var, "_and_", f$var, ".png"
        )
        key <- paste0(y$var, "_by_", x$var, "_and_", f$var)

        all_plots[[key]] <- make_facet_box_plot(
          data        = data,
          x_var       = x$var,
          x_label     = x$label,
          y_var       = y$var,
          y_label     = y$label,
          facet_var   = f$var,
          facet_label = f$label,
          out_dir     = out_dir,
          file_name   = file_name
        )
      }
    }
  }

  return(all_plots)

}
