# Box plots for each (with jittered points)

PlotBoxPlots <- function(data, out_dir, y_vars, x_vars) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  } else {
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Box plot function ---------------------------------------------------------

  make_box_plot <- function(data, x_var, y_var, y_label, group_var, out_dir, file_name) {

    every_nth <- function(n) {
      function(x) x[seq(1, length(x), by = n)]
    }
    
    plot_data <- data %>%
      dplyr::filter(
        !is.na(.data[[x_var]]),
        !is.na(.data[[y_var]])
      ) %>%
      dplyr::mutate(
        !!x_var := if (x_var == "FieldStrength") {
          factor(.data[[x_var]], levels = sort(unique(as.numeric(as.character(.data[[x_var]])))))
        } else {
          factor(
            .data[[x_var]],
            levels = sort(unique(.data[[x_var]]), na.last = TRUE)
          )
        }
      )
    
    # Skip empty plots
    if (nrow(plot_data) == 0) {
      return(NULL)
    }
    
    n_levels <- nlevels(plot_data[[x_var]])
    nth <- max(1, round(n_levels / 50))
    
    # Detect long labels
    label_lengths <- nchar(as.character(plot_data[[x_var]]))
    long_labels <- max(label_lengths, na.rm = TRUE) > 8
    
    # Default logic
    if (n_levels > 10 || long_labels) {
      a <- 45
      j <- 1
    } else {
      a <- 0
      j <- 0.5
    }
    
    # Override specifically for Sequence and ShimMethod
    if (x_var %in% c("Sequence")) {
      a <- 30
      j <- 1
    }
    
    if (x_var %in% c("ShimMethod")) {
      a <- 15
      j <- 1
    }
    
    if (n_levels > 10) {
      plot_width <- 11
    } else if (n_levels > 5) {
      plot_width <- 8
    } else {
      plot_width <- 5
    }
    
    p <- ggplot(
      plot_data,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]]
      )
    ) +
      geom_boxplot(
        fill = "#4C77C2",
        width = 0.65,
        alpha = 0.9,
        color = "black",
        linewidth = 0.4,
        outliers = FALSE
      ) +
      geom_jitter(
        width = 0.06,
        alpha = 0.6,
        color = "black",
        shape = 16,
        size = 0.5
      ) +
      labs(
        title = x_var,
        x = NULL,
        y = bquote(SNR/LW[norm])
      ) +
      scale_x_discrete(
        breaks = every_nth(nth)
      ) +
      theme_comp()

    print(p)

    ggsave(
      filename = file.path(out_dir, file_name),
      plot = p,
      width = plot_width,
      height = 5,
      units = "in",
      dpi = 300,
    )

    return(p)
  }

  # Make box plots ------------------------------------------------------------

  all_boxplots <- list()

  for (y in y_vars) {
    for (x in x_vars) {
      file_name <- paste0("boxplot_", y$var, "_by_", x, ".png")
      key <- paste(y$var, x, sep = "_by_")
      
      all_boxplots[[key]] <- make_box_plot(
        data = data,
        x_var = x,
        y_var = y$var,
        y_label = y$label,
        out_dir = out_dir,
        file_name = file_name
      )
    }
  }
  
  # Remove NULL plots
  all_boxplots <- all_boxplots[!vapply(all_boxplots, is.null, logical(1))]
  
  # One-page combined PDF -----------------------------------------------------
  # A4 width x half A4 height
  page_width  <- 8.27
  page_height <- 11.69 / 2
  
  # Force same font sizes in the combined page
  standardize_plot_for_panel <- function(p) {
    p +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.margin = margin(3, 3, 3, 3),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none"
      )
  }
  
  get_plot <- function(plot_list, name) {
    if (name %in% names(plot_list)) {
      standardize_plot_for_panel(plot_list[[name]])
    } else {
      patchwork::plot_spacer()
    }
  }
  
  for (y in y_vars) {
    
    # Requested order
    p_dp          <- get_plot(all_boxplots, paste(y$var, "DP", sep = "_by_"))
    p_species     <- get_plot(all_boxplots, paste(y$var, "Species", sep = "_by_"))
    
    p_sex         <- get_plot(all_boxplots, paste(y$var, "Sex", sep = "_by_"))
    p_brainregion <- get_plot(all_boxplots, paste(y$var, "VOI", sep = "_by_"))
    p_site        <- get_plot(all_boxplots, paste(y$var, "SiteID", sep = "_by_"))
    
    p_vendor      <- get_plot(all_boxplots, paste(y$var, "Vendor", sep = "_by_"))
    p_field       <- get_plot(all_boxplots, paste(y$var, "FieldStrength", sep = "_by_"))
    p_sequence    <- get_plot(all_boxplots, paste(y$var, "Sequence", sep = "_by_"))
    p_cryoprobe   <- get_plot(all_boxplots, paste(y$var, "Cryoprobe", sep = "_by_"))
    p_shim        <- get_plot(all_boxplots, paste(y$var, "ShimMethod", sep = "_by_"))
    
    keep_y_axis <- function(p) {
      p
    }
    
    remove_y_axis <- function(p) {
      p +
        theme(
          axis.title.y = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank()
        )
    }
    
     row1 <- (p_dp | p_species) +
      patchwork::plot_layout(widths = c(4, 1))
     row1 <- patchwork::wrap_plots(
       keep_y_axis(p_dp),
       remove_y_axis(p_species),
       nrow = 1,
       widths = c(4, 1)
     )
    
    row2 <- (p_sex | p_brainregion | p_site) +
      patchwork::plot_layout(widths = c(0.7, 1.5, 2.5))
    row2 <- patchwork::wrap_plots(
      keep_y_axis(p_sex),
      remove_y_axis(p_brainregion),
      remove_y_axis(p_site),
      nrow = 1,
      widths = c(0.7, 1.5, 2.5)
    )
    
    row3 <- (p_vendor | p_field | p_sequence | p_cryoprobe | p_shim) +
      patchwork::plot_layout(widths = c(2/3, 1.5, 1.5, 2/3, 2/3))
    row3 <- patchwork::wrap_plots(
      keep_y_axis(p_vendor),
      remove_y_axis(p_field),
      remove_y_axis(p_sequence),
      remove_y_axis(p_cryoprobe),
      remove_y_axis(p_shim),
      nrow = 1,
      widths = c(2/3, 1.5, 1.5, 2/3, 2/3)
    )
    
    combined_page <- row1 / row2 / row3 +
      patchwork::plot_layout(heights = c(1, 1, 1))
    combined_page <- patchwork::wrap_plots(
      row1, row2, row3,
      ncol = 1,
      heights = c(1, 1, 1)
    )
    
    base_name <- file.path(out_dir, paste0("combined_boxplots_", y$var, "_onepage"))
    
    # PDF
    ggsave(
      filename = paste0(base_name, ".pdf"),
      plot = combined_page,
      width = page_width,
      height = page_height,
      units = "in",
      dpi = 300,
      device = grDevices::cairo_pdf,
      bg = "transparent"
    )
    
    # PNG (best for slides / presentations)
    ggsave(
      filename = paste0(base_name, ".png"),
      plot = combined_page,
      width = page_width,
      height = page_height,
      units = "in",
      dpi = 300,
      type = "cairo"
    )
    
    # SVG (best for editing in Illustrator / Inkscape)
    ggsave(
      filename = paste0(base_name, ".svg"),
      plot = combined_page,
      width = page_width,
      height = page_height,
      units = "in",
      device = svglite::svglite,
      bg = "transparent"
    )
  }
  
  return(all_boxplots)
}