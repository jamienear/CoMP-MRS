# Pie charts of categorical variables (e.g., species, sex, field strength) to
# show the distribution of DPs across these variables.

PlotPieCharts <- function(data, out_dir, plots_dir, show_pie_charts = TRUE, show_amcharts = TRUE) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  } else {
    unlink(out_dir, recursive = TRUE)
    dir.create(out_dir)
  }
  
  # Custom color palette 
  pie_palette <- c(
    "#4C77C2",  # blue
    "#EA7E2D",  # orange
    "#A7A7A7",  # gray
    "#F2BE00",  # yellow
    "#5B9BD5",  # light blue
    "#70AD47",  # green
    "#264478",  # dark blue
    "#9E480E",  # brown
    "#636363",  # dark gray
    "#997300",  # olive
    "#255E91",  # steel blue
    "#43682B"   # darker green
  )
  
  ### Pie chart function ------------------------------------------------------
  
  make_pie_chart <- function(data, var, plot_title = NULL, file_name = NULL, output_dir) {
    var_sym <- rlang::sym(var)
    
    # AnimalSex should be based on row count, not unique DP-SiteID combinations
    if (var == "AnimalSex") {
      plot_data <- data %>%
        dplyr::filter(!is.na(!!var_sym)) %>%
        dplyr::group_by(!!var_sym) %>%
        dplyr::summarise(n_count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
          perc = n_count / sum(n_count),
          perc_label = scales::percent(perc, accuracy = 0.1),
          legend_label = paste0(as.character(!!var_sym), " (n = ", n_count, ", ", perc_label, ")")
        ) %>%
        dplyr::arrange(dplyr::desc(n_count))
    } else {
      # All other pie charts should be based on unique DP-SiteID combinations
      plot_data <- data %>%
        dplyr::filter(!is.na(!!var_sym), !is.na(DP), !is.na(SiteID)) %>%
        dplyr::distinct(DP, SiteID, !!var_sym) %>%
        dplyr::group_by(!!var_sym) %>%
        dplyr::summarise(n_count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(
          perc = n_count / sum(n_count),
          perc_label = scales::percent(perc, accuracy = 0.1),
          legend_label = paste0(as.character(!!var_sym), " (n = ", n_count, ", ", perc_label, ")")
        ) %>%
        dplyr::arrange(dplyr::desc(n_count))
    }
    
    if (is.null(plot_title)) plot_title <- paste(var, "distribution")
    if (is.null(file_name)) file_name <- paste0(var, "_piechart.png")
    
    fill_vals <- rep(pie_palette, length.out = nrow(plot_data))
    
    p <- ggplot(plot_data, aes(x = 1, y = perc, fill = legend_label)) +
      geom_col(
        color = "white",
        linewidth = 0.4,
        width = 1
      ) +
      geom_text(
        aes(x = 1.3, label = perc_label), # aes(label = count_label),
        position = position_stack(vjust = 0.2),
        color = "white",
        size = 5,
        fontface = "bold",
        family = "serif"
      ) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = fill_vals) +
      xlim(0.5, 1.5) +
      theme_void() +
      labs(
        title = plot_title,
        fill = NULL
      ) +
      theme(
        plot.title = element_text(
          family = "serif",
          face = "bold",
          size = 24,
          hjust = 0.5,
          color = "#4D4D4D"
        ),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text = element_text(
          family = "serif",
          face = "bold",
          size = 12,
          color = "#4D4D4D"
        ),
        legend.key.size = unit(0.35, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        plot.margin = margin(8, 8, 8, 8)
      ) 
    # If the legend has exactly two elements, stack them vertically.
    # Otherwise, keep your original legend arrangement.
    if (nrow(plot_data) <= 4) {
      p <- p +
        guides(
          fill = guide_legend(
            nrow = nrow(plot_data),
            ncol = 1,
            byrow = TRUE
          )
        )
    } else {
      p <- p +
        guides(
          fill = guide_legend(
            nrow = ceiling(nrow(plot_data) / 2),
            byrow = TRUE
          )
        )
    }
    
    
    print(p)
    
    ggsave(
      filename = file.path(output_dir, file_name),
      plot = p,
      height = 5.8,
      width = 5.8,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    
    return(p)
  }
  
  ### Create pie charts with percentages --------------------------------------
  
  pie_specs <- list(
    list(var = "SiteID",         title = "Site (#DPs)"),
    list(var = "Species",  title = "Species (#DPs)"),
    list(var = "Strain",   title = "Animal strain (#DPs)"),
    list(var = "Sex",      title = "Sex (#animals)"),
    list(var = "Vendor",       title = "Vendor (#DPs)"),
    list(var = "FieldStrength",        title = "Field strength (#DPs)"),
    list(var = "Sequence",     title = "Pulse sequence (#DPs)"),
    list(var = "VOI",  title = "Voxel location (#DPs)"),
    list(var = "Averages",     title = "Nr. of averages (#DPs)"),
    list(var = "Cryoprobe",     title = "Cryoprobe (#DPs)"),
    list(var = "ShimMethod",     title = "Shim method (#DPs)")
  )
  
  if (show_pie_charts) {
    
    pie_plots <- lapply(pie_specs, function(x) {
      make_pie_chart(
        data = data,
        var = x$var,
        plot_title = x$title,
        file_name = paste0(x$var, "_piechart.png"),
        output_dir = out_dir
      )
    })
    
    names(pie_plots) <- sapply(pie_specs, `[[`, "var")
    
    ### Combined preview of pie charts ------------------------------------------
    
    combined_pies <- patchwork::wrap_plots(pie_plots, ncol = 4)
    ggsave(
      filename = file.path(plots_dir, "all_piecharts_combined.png"),
      plot = combined_pies,
      width = 20,
      height = 16,
      units = "in",
      dpi = 300
    )
    
  }
  
  ### 3D pie chart function ---------------------------------------------------
  
  # Source - https://stackoverflow.com/a/76142219
  # Posted by Stéphane Laurent
  # Retrieved 2026-03-26, License - CC BY-SA 4.0
  
  make_pie_data_amcharts <- function(data, var) {
    var_sym <- rlang::sym(var)
    
    if (var == "AnimalSex") {
      plot_data <- data %>%
        dplyr::filter(!is.na(!!var_sym)) %>%
        dplyr::group_by(!!var_sym) %>%
        dplyr::summarise(value = dplyr::n(), .groups = "drop")
    } else {
      plot_data <- data %>%
        dplyr::filter(!is.na(!!var_sym), !is.na(DP), !is.na(SiteID)) %>%
        dplyr::distinct(DP, SiteID, !!var_sym) %>%
        dplyr::group_by(!!var_sym) %>%
        dplyr::summarise(value = dplyr::n(), .groups = "drop")
    }
    
    plot_data %>%
      dplyr::rename(group = !!var_sym) %>%
      dplyr::mutate(group = as.character(group))
  }
  
  ### Create 3D pie charts with percentages -----------------------------------
  
  if (show_amcharts) {
    
    amcharts_dir <- file.path(plots_dir, "amcharts_pies")
    dir.create(amcharts_dir, recursive = TRUE, showWarnings = FALSE)
    
    vars <- c(
      "SiteID",
      "Species",
      "Strain",
      "Sex",
      "Vendor",
      "FieldStrength",
      "Sequence",
      "VOI",
      "Averages"
    )
    
    charts <- list()
    
    for (v in vars) {
      dat <- make_pie_data_amcharts(data, v) %>%
        dplyr::arrange(dplyr::desc(value))
      
      charts[[v]] <- amPieChart(
        data = dat,
        category = "group",
        value = "value",
        threeD = TRUE,
        variableDepth = TRUE,
        chartTitle = v,
        legend = TRUE
      )
      
      charts[[v]] <- htmlwidgets::onRender(
        charts[[v]],
        "
    function(el, x) {
      var chart = this.chart;

      if (chart.series && chart.series.values.length > 0) {
        var series = chart.series.values[0];

        series.startAngle = 300;
        series.endAngle   = 660;

        series.slices.template.stroke = am4core.color('#FFFFFF');
        series.slices.template.strokeWidth = 2;
        series.slices.template.strokeOpacity = 1;
      }
    }
    "
      )
      
      print(charts[[v]])
      
      htmlwidgets::saveWidget(
        widget = charts[[v]],
        file = file.path(amcharts_dir, paste0(v, "_piechart.html")),
        selfcontained = TRUE
      )
      
      browseURL(file.path(amcharts_dir,  paste0(v, "_piechart.html")))
    }
    
  }
  
}
