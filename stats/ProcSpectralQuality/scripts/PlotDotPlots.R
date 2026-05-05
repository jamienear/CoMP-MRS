# Dot plots of SNR/LW ratio by subject, before and after outliers removal

PlotDotPlots <- function(data, data_orig, out_dir) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  } else {
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir)
  }
  
  # Dot plot function ---------------------------------------------------------
  
  make_dot_plot <- function(data, x_var, y_var, y_label, group_var, out_dir, file_name) {
    
    every_nth <- function(n) {
      function(x) x[seq(1, length(x), by = n)]
    }
    
    plot_data <- data %>%
      dplyr::filter(
        !is.na(.data[[x_var]]),
        !is.na(.data[[y_var]]),
        !is.na(.data[[group_var]])
      ) %>%
      dplyr::mutate(
        !!x_var := factor(.data[[x_var]], levels = sort(unique(.data[[x_var]])))
      )
    
    p <- ggplot(
      plot_data,
      aes(
        x = .data[[x_var]],
        y = .data[[y_var]]
        # color = factor(.data[[group_var]])
      )
    ) +
      geom_point(
        size = 1,
        color = "#000",
        shape = 16,
        alpha = 0.85
        ) +
      labs(
        title = paste(y_label, "by", x_var),
        x = x_var,
        y = y_label
        # color = group_var
      ) +
      scale_x_discrete(
        breaks = every_nth(12)
      ) +
      theme_comp() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    print(p)
    
    ggsave(
      filename = file.path(out_dir, file_name),
      plot = p,
      width = 7,
      height = 5,
      units = "in",
      dpi = 300
    )
    
    return(p)
  }
  
  # Create dot plots ----------------------------------------------------------
  
  dotplot_orig <- make_dot_plot(
    data = data_orig,
    x_var = "CompID",
    y_var = "SNR_LW_Ratio_norm",
    y_label = "Normalized SNR/LW ratio before outlier removal",
    group_var = "DP",
    out_dir = out_dir,
    file_name = "dotplot_SNR_LW_Ratio_norm_by_Subj_orig.png"
  )
  
  dotplot <- make_dot_plot(
    data = data,
    x_var = "CompID",
    y_var = "SNR_LW_Ratio_norm",
    y_label = "Normalized SNR/LW ratio after outlier removal",
    group_var = "DP",
    out_dir = out_dir,
    file_name = "dotplot_SNR_LW_Ratio_norm_by_Subj.png"
  )
  
}
