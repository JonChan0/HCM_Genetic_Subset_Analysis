#This outlines some general functions used across analyses

log_transformer <- function(data, columns) {
  #Columns are a vector of characters

  if (!identical(columns, '')) {
    for (i in columns) {
      data[[i]] <- log1p(data[[i]])
      data[[i]] <- scale(data[[i]])
      hist(data[[i]], main = i, xlab = i)
    }
  }

  return(data)
}

inverse_rank_transformer <- function(data, columns) {
  #Columns are a vector of characters
  if (!identical(columns, '')) {
    for (i in columns) {
      col_label <- i
      col_label_before <- columns[match(i, columns) - 1]
      mask <- !is.na(data[[i]]) #Mask to filter out NA values

      # Filter out NA values for that phenotype
      data_copy <- subset(data, mask)
      data_copy[[i]] <- RNOmni::RankNorm(data_copy[[i]]) #Performing rank-inverse normalisation
      data_copy <- select(data_copy, 'HCR_IDs', all_of(i)) #Taking only the column with HCR_ID and the column of interest

      data <- select(data, -i) %>% #Removing the original column and adding back the new column after rank-inverse normalisation
        left_join(data_copy, by = 'HCR_IDs') %>%
        relocate(i, .after = !!col_label_before) #Moving the appended column back to its original index

      hist(data[[i]], main = i, xlab = i)
    }
  }

  return(data)
}

#--------------------------------------------------------------------------------
# This function performs group comparison plots + statistical tests for an input_param_vector of c(var, label).

plot_tester <- function(
  input_tb,
  input_param_vector,
  group1_name = NULL,
  group2_name = NULL,
  output_path,
  continuous_or_discrete = "continuous",
  groups = NULL, # Groups refers to the c(label1, label2, label3...)
  ref_group = NULL,
  min_cell_n = 10
) {
  phenotype_colname <- input_param_vector[[1]]
  xlabel <- input_param_vector[[2]]

  if (is.null(groups)) {
    if (!is.null(group1_name) && !is.null(group2_name)) {
      groups <- c(group1_name, group2_name)
    } else {
      groups <- sort(unique(input_tb$group))
    }
  }
  groups <- unique(groups)

  plot_tb <- input_tb %>%
    filter(!is.na(.data$group)) %>%
    filter(.data$group %in% groups) %>%
    filter(!is.na(.data[[phenotype_colname]])) %>%
    mutate(group = factor(.data$group, levels = groups))

  group_counts <- plot_tb %>%
    count(.data$group, name = "n")
  group_count_label <- group_counts %>%
    mutate(label = str_c(as.character(.data$group), " (n=", .data$n, ")")) %>%
    pull(.data$label) %>%
    str_c(collapse = ", ")

  if (continuous_or_discrete == "continuous") {
    if (length(groups) < 2) {
      stop("Need at least 2 groups to compare.")
    }

    if (length(groups) == 2) {
      vals_g1 <- filter(plot_tb, .data$group == groups[[1]])
      vals_g2 <- filter(plot_tb, .data$group == groups[[2]])

      wilcox_test_results <- wilcox.test(
        vals_g1[[phenotype_colname]],
        vals_g2[[phenotype_colname]],
        conf.int = TRUE
      )

      med1 <- stats::median(vals_g1[[phenotype_colname]], na.rm = TRUE)
      med2 <- stats::median(vals_g2[[phenotype_colname]], na.rm = TRUE)

      # Report group1 - group2 (not group2 - group1)
      diff_median <- med1 - med2

      test_caption <- str_c(
        "Nominal 2-sided Wilcoxon p-value = ",
        as.character(signif(wilcox_test_results$p.value, digits = 3)),
        "\nMedian(",
        groups[[1]],
        ") - Median(",
        groups[[2]],
        ") = ",
        as.character(signif(diff_median, digits = 3))
      )
    } else {
      kw <- kruskal.test(
        stats::as.formula(str_c(phenotype_colname, " ~ group")),
        data = plot_tb
      )
      test_caption <- str_c(
        "Nominal Kruskal-Wallis p-value = ",
        as.character(signif(kw$p.value, digits = 3))
      )
    }

    dens_plot <- ggplot(plot_tb) +
      geom_density(
        aes(x = .data[[phenotype_colname]], fill = .data$group),
        alpha = 0.3
      ) +
      xlab(xlabel) +
      ylab("Density") +
      labs(
        title = str_wrap(
          str_c(
            "Distribution of",
            xlabel,
            "across groups:",
            group_count_label,
            sep = " "
          ),
          width = 80
        ),
        caption = test_caption,
        fill = "Group"
      )
    print(dens_plot)
    ggsave(
      filename = str_c(
        output_path,
        "/",
        str_c(phenotype_colname, "density_plot.png", sep = "_")
      ),
      plot = dens_plot,
      dpi = 600
    )

    hist_plot <- ggplot(plot_tb) +
      geom_histogram(
        aes(x = .data[[phenotype_colname]], fill = .data$group),
        alpha = 0.5
      ) +
      xlab(xlabel) +
      ylab("Count") +
      labs(
        title = str_wrap(
          str_c(
            "Histogram of",
            xlabel,
            "across groups:",
            group_count_label,
            sep = " "
          ),
          width = 80
        ),
        caption = test_caption,
        fill = "Group"
      )
    print(hist_plot)
    ggsave(
      filename = str_c(
        output_path,
        "/",
        str_c(phenotype_colname, "hist_plot.png", sep = "_")
      ),
      plot = hist_plot,
      dpi = 600
    )

    if (length(groups) == 2) {
      output_test <- tribble(
        ~pheno                       , ~group1                         , ~group2                         , ~pval                       , ~estimate , ~lowerCI , ~upperCI ,
        ~median_group1               , ~median_group2                  , ~diff_median                    ,
        phenotype_colname            , groups[[1]]                     , groups[[2]]                     , wilcox_test_results$p.value ,
        wilcox_test_results$estimate , wilcox_test_results$conf.int[1] , wilcox_test_results$conf.int[2] ,
        med1                         , med2                            , diff_median
      )
    } else {
      output_test <- tribble(
        ~pheno            , ~test            , ~pval      ,
        phenotype_colname , "Kruskal-Wallis" , kw$p.value
      )

      if (!is.null(ref_group)) {
        if (!ref_group %in% groups) {
          stop("ref_group must be one of `groups`.")
        }
        comps <- setdiff(groups, ref_group)
        pairwise_tb <- purrr::map_dfr(comps, function(g) {
          v_ref <- filter(plot_tb, .data$group == ref_group)
          v_g <- filter(plot_tb, .data$group == g)
          wt <- wilcox.test(
            v_ref[[phenotype_colname]],
            v_g[[phenotype_colname]],
            conf.int = TRUE
          )

          med_ref <- stats::median(v_ref[[phenotype_colname]], na.rm = TRUE)
          med_g <- stats::median(v_g[[phenotype_colname]], na.rm = TRUE)

          # IMPORTANT: group1 - group2, where group1 is ref_group and group2 is g
          diff_median <- med_ref - med_g

          tibble(
            pheno = phenotype_colname,
            group1 = ref_group,
            group2 = g,
            pval = wt$p.value,
            estimate = wt$estimate,
            lowerCI = wt$conf.int[1],
            upperCI = wt$conf.int[2],
            median_group1 = med_ref,
            median_group2 = med_g,
            diff_median = diff_median
          )
        })
        output_test <- list(global = output_test, pairwise_vs_ref = pairwise_tb)
      }
    }

    return(output_test)
  }

  if (continuous_or_discrete == "discrete") {
    if (length(groups) < 2) {
      stop("Need at least 2 groups to compare.")
    }

    contingency_table <- table(plot_tb$group, plot_tb[[phenotype_colname]])

    keep_cols <- apply(contingency_table, 2, function(x) all(x > min_cell_n))
    filtered_table <- contingency_table[, keep_cols, drop = FALSE]

    if (ncol(filtered_table) < 2) {
      stop("After filtering by min_cell_n, <2 phenotype levels remain.")
    }

    chisq_results <- chisq.test(filtered_table)

    heatmap_plot <- ggplot(
      plot_tb,
      aes(x = .data$group, y = .data[[phenotype_colname]])
    ) +
      geom_bin2d() +
      stat_bin2d(geom = "text", aes(label = ..count..), col = "white") +
      scale_fill_viridis(name = "Count") +
      xlab("Group") +
      ylab(xlabel) +
      labs(
        title = str_wrap(
          str_c(
            "Heatmap of",
            xlabel,
            "across groups:",
            group_count_label,
            sep = " "
          ),
          width = 80
        ),
        caption = str_c(
          "Nominal Chi-squared p-value (post-filtering for cell sizes > ",
          min_cell_n,
          ") = ",
          as.character(signif(chisq_results$p.value, digits = 3))
        )
      )

    print(heatmap_plot)
    ggsave(
      filename = str_c(
        output_path,
        "/",
        str_c(phenotype_colname, "heatmap_plot.png", sep = "_")
      ),
      plot = heatmap_plot,
      dpi = 600,
      width = 10,
      height = 6
    )

    output_chisq_test <- tribble(
      ~pheno            , ~test         , ~pval                 ,
      phenotype_colname , "Chi-squared" , chisq_results$p.value
    )

    return(output_chisq_test)
  }

  stop("continuous_or_discrete must be either 'continuous' or 'discrete'.")
}

#This provides a summary plot of many statistical tests in a single summary scatter plot
#This is for the rare variant subset analysis comparing across different discrete groups
summary_plotter <- function(
  input_tb,
  output_path = "../output/plots/1_MAIN/summary_plots/1_MAIN/",
  suffix = '',
  mtc = "bonferroni",
  pval_col = 'pval',
  confint = FALSE,
  se = FALSE,
  label = FALSE,
  poster = FALSE,
  test = "Wilcoxon",
  shape_values = "",
  xlabel = "Estimate",
  width = 9,
  height = 6,
  comporder = c(
    "Thick vs. Thin",
    "MYBPC3 vs. MYH7",
    "MYBPC3_pLOF vs. MYBPC3_nonLOF",
    "P/LP vs. VUS"
  )
) {
  # ---- Standardise expected columns across "group comparisons" and "continuous predictor" results ----
  # group comparisons typically have: comparison, group1, group2
  # continuous predictor tests typically have: predictor (and no comparison)
  is_group_comp <- "comparison" %in% names(input_tb)

  # Safe column accessor: returns an all-NA vector if the column is absent
  col_or_na <- function(.df, .name) {
    if (.name %in% names(.df)) .df[[.name]] else rep(NA, nrow(.df))
  }

  plot_tb <- input_tb %>%
    mutate(
      pheno = dplyr::coalesce(col_or_na(., "pheno"), col_or_na(., "variable")),
      plot_type = ifelse(
        is_group_comp,
        "group_comparison",
        "continuous_predictor"
      ),
      Label = dplyr::coalesce(col_or_na(., "Label")),
      estimate = dplyr::coalesce(
        col_or_na(., "estimate"),
        col_or_na(., "Estimate")
      ),
      pval = dplyr::coalesce(
        col_or_na(., pval_col)
      ),
      lowerCI = dplyr::coalesce(
        col_or_na(., "lowerCI"),
        col_or_na(., "lower_ci"),
        col_or_na(., "conf.low"),
        col_or_na(., "estimate") - 1.96 * col_or_na(., "SE")
      ),
      upperCI = dplyr::coalesce(
        col_or_na(., "upperCI"),
        col_or_na(., "upper_ci"),
        col_or_na(., "conf.high"),
        col_or_na(., "estimate") + 1.96 * col_or_na(., "SE")
      ),
      lowerSE = dplyr::coalesce(col_or_na(., "lowerSE")),
      upperSE = dplyr::coalesce(col_or_na(., "upperSE"))
    ) %>%
    filter(!is.na(.data$pheno), !is.na(.data$estimate), !is.na(.data$pval))

  if (isFALSE(is_group_comp)) {
    plot_tb <- plot_tb %>%
      mutate(
        comparison = ifelse(
          "predictor" %in% names(.),
          as.character(.data$predictor),
          "predictor"
        )
      )
  }

  if (is_group_comp) {
    plot_tb <- plot_tb %>%
      mutate(comparison = factor(.data$comparison, levels = comporder))
  } else {
    plot_tb <- plot_tb %>%
      mutate(
        comparison = factor(.data$comparison, levels = unique(.data$comparison))
      )
  }

  # ---- Multiple testing correction + labels ----
  using_adj_p <- identical(pval_col, "adj_p")

  plot_tb <- plot_tb %>%
    mutate(
      pval = dplyr::case_when(
        using_adj_p ~ .data$pval, # If using adj_p, don't do any further adjustment; treat .data$pval as already adjusted
        mtc == "fdr" ~ p.adjust(.data$pval, "fdr"),
        mtc == "BY" ~ p.adjust(.data$pval, "BY"),
        mtc == "holm" ~ p.adjust(.data$pval, "holm"),
        TRUE ~ .data$pval
      )
    ) %>%
    mutate(
      sig = dplyr::case_when(
        using_adj_p ~ .data$pval < 0.05,
        mtc == "bonferroni" ~ .data$pval < (0.05 / nrow(.)),
        mtc %in% c("fdr", "BY", "holm") ~ .data$pval < 0.05,
        TRUE ~ .data$pval < 0.05
      )
    ) %>%
    mutate(
      sig_label = dplyr::case_when(
        !isTRUE(label) ~ "",
        .data$sig & 'Label_short' %in% colnames(plot_tb) ~ .data$Label_short,
        .data$sig ~ .data$Label,
        TRUE ~ ""
      )
    )

  View(plot_tb)

  ylab_text <- if (using_adj_p || mtc %in% c("fdr", "BY", "holm")) {
    "-log10(Adjusted p-value)"
  } else {
    "-log10(p-value)"
  }

  # ---- Plot ----
  summary_cont_plot <- ggplot(
    plot_tb,
    aes(
      x = .data$estimate,
      y = -log10(.data$pval),
      col = .data$Label,
      shape = .data$comparison
    )
  ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_point() +
    {
      if ("colour" %in% names(plot_tb)) {
        colour_map <- plot_tb |>
          dplyr::select(.data$Label, .data$colour) |>
          dplyr::distinct() |>
          tidyr::drop_na()

        ggplot2::scale_colour_manual(
          values = rlang::set_names(colour_map$colour, colour_map$Label),
          labels = scales::label_wrap(30)
        )
      }
    } +
    labs(
      shape = ifelse(is_group_comp, "Comparison", "Predictor"),
      col = "Phenotype"
    ) +
    scale_x_continuous(n.breaks = 10) +
    xlab(xlabel) +
    ylab(ylab_text) +
    theme(legend.key.spacing.y = grid::unit(8, "points")) +
    labs(
      title = stringr::str_wrap(
        stringr::str_c(
          "Summary plot of ",
          test,
          " results for ",
          length(unique(plot_tb$pheno)),
          " phenotypes across ",
          length(unique(plot_tb$comparison)),
          ifelse(is_group_comp, " group comparisons", " predictors")
        )
      )
    )

  if (!identical(shape_values, "")) {
    summary_cont_plot <- summary_cont_plot +
      scale_shape_manual(values = shape_values)
  }

  if (isTRUE(label)) {
    summary_cont_plot <- summary_cont_plot +
      ggrepel::geom_text_repel(aes(label = .data$sig_label), force = 5)
  }

  if (isFALSE(is_group_comp)) {
    summary_cont_plot <- summary_cont_plot +
      # Remove only the Predictor legend sub-box, not the Pheno one
      guides(shape = "none")
  }

  if (using_adj_p) {
    # No additional lines/captions from adjustment method; already adjusted p-values supplied
    summary_cont_plot <- summary_cont_plot +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      labs(caption = "Adjusted p-values provided in `adj_p`")
  } else if (mtc == "bonferroni") {
    summary_cont_plot <- summary_cont_plot +
      geom_hline(
        yintercept = -log10(0.05 / nrow(plot_tb)),
        linetype = "dashed"
      ) +
      labs(
        caption = stringr::str_c(
          "Bonferroni threshold = ",
          signif(0.05 / nrow(plot_tb), 3)
        )
      )
  } else if (mtc %in% c("fdr", "BY", "holm")) {
    summary_cont_plot <- summary_cont_plot +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      labs(
        caption = stringr::str_c(
          ifelse(mtc == "fdr", "FDR BH", pval),
          " p-value correction applied"
        )
      )
  }

  append <- suffix

  if (isTRUE(confint)) {
    summary_cont_plot <- summary_cont_plot +
      geom_errorbar(
        aes(xmin = .data$lowerCI, xmax = .data$upperCI),
        alpha = 0.5,
        width = 0.05
      )
    append <- str_c(suffix, "_confint")
  } else if (isTRUE(se)) {
    summary_cont_plot <- summary_cont_plot +
      geom_errorbar(
        aes(xmin = .data$lowerSE, xmax = .data$upperSE),
        alpha = 0.5,
        width = 0.05
      )
    append <- str_c(suffix, "_SE")
  }

  if (isTRUE(poster)) {
    summary_cont_plot <- summary_cont_plot +
      theme(
        legend.position = "bottom",
        legend.box = "vertical",
        text = element_text(size = 16)
      ) +
      scale_colour_brewer(palette = "Dark2", label = scales::label_wrap(20)) +
      geom_point(size = 5) +
      labs(title = "", caption = "")

    if (isTRUE(confint)) {
      summary_cont_plot <- summary_cont_plot +
        geom_errorbar(
          aes(xmin = .data$lowerCI, xmax = .data$upperCI),
          linewidth = 1,
          alpha = 0.75
        )
    } else if (isTRUE(se)) {
      summary_cont_plot <- summary_cont_plot +
        geom_errorbar(
          aes(xmin = .data$lowerSE, xmax = .data$upperSE),
          linewidth = 1,
          alpha = 0.75
        )
    }

    append <- stringr::str_c(append, "_poster")
  }

  summary_cont_plot %>% print()

  ggsave(
    filename = stringr::str_c(
      output_path,
      test,
      "_summary_plot_",
      if (using_adj_p) "adj_p" else mtc,
      append,
      ".png"
    ),
    plot = summary_cont_plot,
    dpi = 600,
    width = width,
    height = if (isTRUE(poster)) 12 else height
  )

  ggsave(
    filename = stringr::str_c(
      output_path,
      test,
      "_summary_plot_",
      if (using_adj_p) "adj_p" else mtc,
      append,
      ".tiff"
    ),
    plot = summary_cont_plot,
    dpi = 600,
    width = width,
    height = if (isTRUE(poster)) 12 else height
  )

  invisible(summary_cont_plot)
}

#This plots a forest plot of significant associations from the linear regression results
summary_forest_plotter <- function(
  #N.B Label is used as the facet variable btw! and predictor_label is the y-axis labels
  lm_results_scaled,
  output_dir,
  facet_order,
  predictor_labels = c(cv_prs = "PRS", sarcomere = "Sarcomere Status"),
  alpha = 0.05,
  file_name = "forest_plot.png",
  width = 12,
  height = 6,
  dpi = 600,
  sigphenos_only = T
) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Identify significant associations (pheno x predictor) based on adjusted p
  sig_phenos <- lm_results_scaled |>
    filter(adj_p < alpha) |>
    select(pheno, predictor_reported) |>
    mutate(combo = str_c(pheno, predictor_reported, sep = "_"))

  if (nrow(sig_phenos) == 0) {
    stop("No significant associations found at alpha = ", alpha, ".")
  }

  # Forest plot table (scaled betas)
  forest_tb <- lm_results_scaled |>
    mutate(combo = str_c(pheno, predictor_reported, sep = "_"))

  if (sigphenos_only == TRUE) {
    forest_tb <- forest_tb |>
      filter(combo %in% sig_phenos$combo)
  }

  forest_tb <- forest_tb |>
    mutate(
      conf_low = estimate - 1.96 * SE,
      conf_high = estimate + 1.96 * SE,
      sig_label = case_when(
        adj_p <= 0.001 ~ "***",
        adj_p <= 0.01 ~ "**",
        adj_p <= 0.05 ~ "*",
        TRUE ~ ""
      ),
      # Maps each value in `predictor_reported` to a human-readable label using `predictor_labels`;
      # if no mapping exists for a value, it keeps the original value unchanged (`.default = predictor_reported`).
      predictor_label = dplyr::recode(
        predictor_reported,
        !!!predictor_labels,
        .default = predictor_reported
      ),
      predictor_label = factor(predictor_label, levels = facet_order)
    )

  # ---- Alignment: enforce identical y ordering within each facet ----
  y_levels_by_facet <- forest_tb |>
    distinct(predictor_label, Label) |>
    group_by(predictor_label) |>
    mutate(.y_id = row_number()) |>
    ungroup()

  forest_tb <- forest_tb |>
    left_join(y_levels_by_facet, by = c("predictor_label", "Label")) |>
    mutate(Label_facet = fct_reorder(Label, .y_id))

  p_forest <- ggplot(
    forest_tb,
    aes(x = estimate, y = Label_facet, fill = Label)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.2) +
    scale_fill_manual(
      values = setNames(forest_tb$colour, forest_tb$Label),
      guide = "none"
    ) +
    geom_point(shape = 23, size = 2) +
    geom_text(aes(label = sig_label, x = estimate), vjust = -0.25, size = 5) +
    facet_wrap(~predictor_label, ncol = 1) +
    labs(
      x = "Phenotype Difference (SD)",
      y = NULL
    ) +
    theme(strip.text = element_text(face = "bold"))

  ggsave(
    filename = file.path(output_dir, file_name),
    plot = p_forest,
    width = width,
    height = height,
    dpi = dpi
  )

  p_forest
}

#---------------------------------------------------------------------------------------
#This function performs association analysis using linear/logistic regression to enable adjustment for covariates
linear_logistic_regression_modeller <- function(
  input_tb,
  input_param_vector,
  group1_name = NULL,
  group2_name = NULL,
  covars = c("age", "sex"),
  continuous_or_discrete = "continuous",
  model_return = FALSE,
  summary_print = FALSE,
  predictor = "group",
  predictor_type = c("group", "continuous"),
  additional_predictors_to_report = character(),
  bsa_phenos = c(
    "wallthkmax",
    "echomax",
    "rvs_amp_avg_glob",
    "rvs_rate_dia_avg_glob",
    "rvs_rate_sys_avg_glob",
    "rvs_ttp_avg_glob"
  )
) {
  predictor_type <- match.arg(predictor_type)

  pheno_of_interest <- input_param_vector[[1]]
  plot_tb <- input_tb

  # ---- Validate columns ----
  required_vars <- unique(c(
    pheno_of_interest,
    covars,
    predictor,
    additional_predictors_to_report
  ))
  missing_vars <- setdiff(required_vars, names(plot_tb))
  if (length(missing_vars) > 0) {
    stop(
      "Missing column(s) in input_tb: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # ---- Predictor handling (main predictor) ----
  # If predictor is a group-like categorical variable:
  # - if group1_name & group2_name are provided: enforce 2-level comparison
  # - otherwise: use all observed levels of the predictor (drop NA)
  if (predictor_type == "group") {
    has_group_pair <- !is.null(group1_name) && !is.null(group2_name)

    if (has_group_pair) {
      plot_tb <- plot_tb |>
        dplyr::filter(.data[[predictor]] %in% c(group1_name, group2_name)) |>
        dplyr::mutate(
          !!predictor := factor(
            .data[[predictor]],
            levels = c(group1_name, group2_name)
          )
        )
    } else {
      plot_tb <- plot_tb |>
        dplyr::filter(!is.na(.data[[predictor]])) |>
        dplyr::mutate(!!predictor := factor(.data[[predictor]]))
    }
  }

  # (Optional) coerce additional predictors that are character to factors.
  # This supports "either continuous or group" without requiring a separate type map.
  if (length(additional_predictors_to_report) > 0) {
    plot_tb <- plot_tb |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(additional_predictors_to_report),
          ~ if (is.character(.x)) factor(.x) else .x
        )
      )
  }

  # ---- Build model ----
  include_bsa <- pheno_of_interest %in% bsa_phenos
  rhs_terms <- c(covars, if (include_bsa) "bsa", predictor)

  model_formula <- stats::as.formula(
    stringr::str_c(
      pheno_of_interest,
      " ~ ",
      stringr::str_c(rhs_terms, collapse = " + ")
    )
  )

  print(model_formula)

  model <- if (continuous_or_discrete == "continuous") {
    stats::lm(model_formula, data = plot_tb, na.action = stats::na.omit)
  } else if (continuous_or_discrete == "discrete") {
    stats::glm(
      model_formula,
      family = stats::binomial(),
      data = plot_tb,
      na.action = stats::na.omit
    )
  } else {
    stop("continuous_or_discrete must be either 'continuous' or 'discrete'.")
  }

  if (isTRUE(summary_print)) {
    summary(model)
  }

  # ---- Extract coefficients to report ----
  # Include the main predictor plus any requested additional predictors.
  predictors_to_report <- unique(c(predictor, additional_predictors_to_report))

  # Helper: which terms correspond to a predictor?
  # - for numeric predictors it's an exact match (e.g. "age")
  # - for factors it's multiple rows like "groupB" / "sexMale"
  match_predictor_terms <- function(.terms, .pred) {
    .terms == .pred | stringr::str_starts(.terms, stringr::str_c(.pred))
  }

  tidy_tb <- broom::tidy(model)

  coef_tb <- purrr::map_dfr(predictors_to_report, function(pred) {
    pred_terms <- tidy_tb |>
      dplyr::filter(match_predictor_terms(.data$term, pred)) |>
      dplyr::transmute(
        predictor_reported = pred,
        term = .data$term,
        estimate = .data$estimate,
        SE = .data$std.error,
        pval = .data$p.value
      )

    pred_terms
  })

  if (nrow(coef_tb) < 1) {
    stop(
      "Could not identify any requested predictor term(s) in the model. Check names in `predictor` / `additional_predictors_to_report`."
    )
  }

  output_tb <- tibble::tibble(
    pheno = pheno_of_interest,
    predictor = predictor,
    predictor_type = predictor_type,
    group1 = if (predictor_type == "group") group1_name else NA_character_,
    group2 = if (predictor_type == "group") group2_name else NA_character_
  ) |>
    dplyr::slice(rep(1, nrow(coef_tb))) |>
    dplyr::bind_cols(coef_tb)

  if (isTRUE(model_return)) model else output_tb
}

#This plots diagnostic plots of the linear/logistic regression model to check assumptions
diagnostic_plotter <- function(input_model, pheno) {
  broom::augment(input_model)

  par(mfrow = c(2, 2))
  plot(input_model, main = pheno)
}

#This hardcodes a colour palette into the output tibble from linear regression.

colour_hardcoder <- function(input_tb) {
  #Now assign a distinct colour hex code for each pheno in the output_tb, taking the HEX code from the colour palette of Set2 from ColourBrewer if <8 phenotypes, and from Paired if >=8 phenotypes
  # Assign colours by phenotype in a stable way:
  # - build the mapping on ALL phenotypes you want consistently coloured (not just those present in this output_tb)
  # - then join the colour back in

  phenos_all <- unique(input_tb$pheno)
  n_pheno <- length(phenos_all)

  if (n_pheno <= 8) {
    palette <- RColorBrewer::brewer.pal(8, "Set2")[seq_len(n_pheno)]
  } else if (n_pheno <= 12) {
    palette <- RColorBrewer::brewer.pal(12, "Paired")[seq_len(n_pheno)]
  } else {
    palette <- grDevices::hcl.colors(n_pheno, palette = "Dynamic")
  }

  pheno_colour_mapping <- tibble::tibble(
    pheno = phenos_all,
    colour = palette
  )

  output_tb <- input_tb |>
    dplyr::left_join(pheno_colour_mapping, by = "pheno")

  return(output_tb)
}


#-------------------------------------------------------------
#Rare Variant Subset Analysis x Endophenotype Functions

#This function generates each tibble for each comparison
comparison_pheno_tb_gen <- function(
  phenotype_tb,
  param_list,
  rarevar_class_tb,
  var_class = 'primary',
  only_singles = T
) {
  #This defines the comparison groups based on the param_list input e.g. for thick vs. thin, G1 = Thick, G2 = Thin
  group1_classes <- param_list[[1]]
  group1_name <- param_list[[2]]
  group2_classes <- param_list[[3]]
  group2_name <- param_list[[4]]

  if (var_class == 'primary') {
    # When using primary classification based on gene-based

    output_tb <- phenotype_tb %>%
      mutate(
        group = case_when(
          class_primary %in% group1_classes ~ group1_name,
          class_primary %in% group2_classes ~ group2_name,
          T ~ NA
        )
      )
  } else if (var_class == 'secondary') {
    # When using secondary classification based on pathogenicity

    output_tb <- phenotype_tb %>%
      mutate(
        group = case_when(
          class_secondary %in% group1_classes ~ group1_name,
          class_secondary %in% group2_classes ~ group2_name,
          T ~ NA
        )
      )
  }

  output_tb <- mutate(
    output_tb,
    group = factor(group, levels = c(group2_name, group1_name))
  ) %>%
    filter(!is.na(group)) #Filter to only individuals in either group

  #This filters down to single dosage individuals,
  if (isTRUE(only_singles)) {
    if (var_class == 'primary') {
      single_dosage_individuals <- rarevar_class_tb %>%
        mutate(
          total_dosage = select(., MYBPC3_pLOF:FHOD3) %>% rowSums(na.rm = T)
        ) %>%
        filter(total_dosage == 1)
    } else if (var_class == 'secondary') {
      single_dosage_individuals <- rarevar_class_tb %>%
        mutate(
          total_dosage = select(., MYBPC3_P_LP:FHOD3_VUS_FP_VUS) %>%
            rowSums(na.rm = T)
        ) %>%
        filter(total_dosage == 1)
    }
    output_tb <- filter(
      output_tb,
      HCR_IDs %in% single_dosage_individuals$HCR_ID
    )
  }
  return(output_tb)
}

#------------------------------------------------------------------------------------------------
#For meta-analysis
hcmr_significant_results_grabber <- function(
  input_comparison_tb,
  group1,
  group2,
  sig_phenotypes,
  mean = T
) {
  if (isTRUE(mean)) {
    output_tb <- input_comparison_tb %>%
      summarise(
        Study = 'HCMR',
        Group1 = group1,
        Group2 = group2,
        Phenotype = sig_phenotypes,
        N1 = map_int(
          sig_phenotypes,
          ~ sum(!is.na(filter(input_comparison_tb, group == group1)[[.]]))
        ),
        N2 = map_int(
          sig_phenotypes,
          ~ sum(!is.na(filter(input_comparison_tb, group == group2)[[.]]))
        ),
        Mean1 = map_dbl(
          sig_phenotypes,
          ~ mean(filter(input_comparison_tb, group == group1)[[.]], na.rm = T)
        ),
        Mean2 = map_dbl(
          sig_phenotypes,
          ~ mean(filter(input_comparison_tb, group == group2)[[.]], na.rm = T)
        ),
        SD1 = map_dbl(
          sig_phenotypes,
          ~ sd(filter(input_comparison_tb, group == group1)[[.]], na.rm = T)
        ),
        SD2 = map_dbl(
          sig_phenotypes,
          ~ sd(filter(input_comparison_tb, group == group2)[[.]], na.rm = T)
        )
      )
  } else {
    output_tb <- input_comparison_tb %>%
      summarise(
        Study = 'HCMR',
        Group1 = group1,
        Group2 = group2,
        Phenotype = sig_phenotypes,
        N1 = map_int(
          sig_phenotypes,
          ~ sum(!is.na(filter(input_comparison_tb, group == group1)[[.]]))
        ),
        N2 = map_int(
          sig_phenotypes,
          ~ sum(!is.na(filter(input_comparison_tb, group == group2)[[.]]))
        ),
        Median1 = map_dbl(
          sig_phenotypes,
          ~ median(filter(input_comparison_tb, group == group1)[[.]], na.rm = T)
        ),
        Median2 = map_dbl(
          sig_phenotypes,
          ~ median(filter(input_comparison_tb, group == group2)[[.]], na.rm = T)
        ),
        FQ_1 = map_dbl(
          sig_phenotypes,
          ~ quantile(
            filter(input_comparison_tb, group == group1)[[.]],
            0.25,
            na.rm = T
          )
        ),
        TQ_1 = map_dbl(
          sig_phenotypes,
          ~ quantile(
            filter(input_comparison_tb, group == group1)[[.]],
            0.75,
            na.rm = T
          )
        ),
        FQ_2 = map_dbl(
          sig_phenotypes,
          ~ quantile(
            filter(input_comparison_tb, group == group2)[[.]],
            0.25,
            na.rm = T
          )
        ),
        TQ_2 = map_dbl(
          sig_phenotypes,
          ~ quantile(
            filter(input_comparison_tb, group == group2)[[.]],
            0.75,
            na.rm = T
          )
        )
      )
  }
  return(output_tb)
}

meta_cont_from_tb <- function(
  meta_tb,
  pheno,
  group1,
  group2,
  plot_output_path = 'genetic_subset_analysis/RareVar/3_output/plots/1_Endophenotype/3_Meta/MD/CMR/',
  excluded_studies = '',
  effectsize = 'MD',
  tau_method = 'REML',
  remove_medians = F,
  method_to_estimate_mean = 'Luo',
  prediction = T,
  random = F
) {
  filtered_tb <- meta_tb %>%
    filter(Phenotype == pheno, Group1 == group1, Group2 == group2)

  if (isFALSE(remove_medians)) {
    meta_cont_result <- metacont(
      n.e = filtered_tb$N1,
      n.c = filtered_tb$N2,
      mean.e = filtered_tb$Mean1,
      mean.c = filtered_tb$Mean2,
      sd.e = filtered_tb$SD1,
      sd.c = filtered_tb$SD2,
      median.e = filtered_tb$Median1,
      median.c = filtered_tb$Median2,
      q1.e = filtered_tb$FQ_1,
      q3.e = filtered_tb$TQ_1,
      q1.c = filtered_tb$FQ_2,
      q3.c = filtered_tb$TQ_2,
      method.mean = method_to_estimate_mean,
      studlab = filtered_tb$Study,
      exclude = filtered_tb$Study %in% excluded_studies,
      sm = effectsize,
      method.smd = 'Hedges',
      common = T,
      random = random,
      prediction = T,
      method.tau = tau_method,
      title = str_c(
        'Meta-analysis of',
        nrow(filtered_tb),
        'studies comparing',
        pheno,
        'in',
        group1,
        'vs.',
        group2,
        sep = ' '
      ),
      label.e = group1,
      label.c = group2
    )
  } else {
    filtered_tb <- filter(filtered_tb, is.na(Median1))

    meta_cont_result <- metacont(
      n.e = filtered_tb$N1,
      n.c = filtered_tb$N2,
      mean.e = filtered_tb$Mean1,
      mean.c = filtered_tb$Mean2,
      sd.e = filtered_tb$SD1,
      sd.c = filtered_tb$SD2,
      studlab = filtered_tb$Study,
      exclude = filtered_tb$Study %in% excluded_studies,
      sm = effectsize,
      method.smd = 'Hedges',
      common = T,
      random = random,
      prediction = T,
      method.tau = tau_method,
      title = str_c(
        'Meta-analysis of',
        nrow(filtered_tb),
        'studies comparing',
        pheno,
        'in',
        group1,
        'vs.',
        group2,
        sep = ' '
      ),
      label.e = group1,
      label.c = group2
    )
  }

  print(summary(meta_cont_result))

  #Output a Forest plot
  output_name <- str_c(
    plot_output_path,
    '/forest/',
    'meta_',
    group1,
    '_',
    group2,
    pheno,
    '_',
    '_forest.png'
  )
  png(filename = output_name, width = 7200, height = 3600, res = 600)

  forest_plot <- forest(
    meta_cont_result,
    common = T,
    random = random,
    overall = T,
    prediction = F,
    print.Q = F,
    print.pval.Q = T,
    test.overall = T
  )
  dev.off()

  #Output funnel plot
  output_name2 <- str_c(
    plot_output_path,
    '/funnel/',
    'meta_',
    group1,
    '_',
    group2,
    pheno,
    '_funnel.png'
  )
  png(filename = output_name2, width = 7200, height = 2600, res = 600)
  funnel_plot <- funnel(meta_cont_result, studlab = T)
  dev.off()

  #Print plots
  forest_plot <- forest(
    meta_cont_result,
    common = T,
    random = random,
    overall = T,
    prediction = F,
    print.Q = F,
    print.pval.Q = T,
    test.overall = T
  )

  funnel_plot <- funnel(meta_cont_result, studlab = T)

  return(meta_cont_result)
}

#------------------------------------------------------------------------------------------------
#Publication summary figures for the rare variant subset comparisons

#This plots a panelled summary scatter (volcano) plot with one panel per genotypic comparison,
#sharing the same x-axis (coefficient estimate) and y-axis (-log10 adjusted p-value) across all panels
summary_panel_plotter <- function(
  input_tb,
  output_path,
  file_name = 'LinearRegression_summary_panelplot_fdr.png',
  comporder = c(
    'Thick vs. Thin',
    'MYBPC3 vs. MYH7',
    'MYBPC3_pLOF vs. MYBPC3_nonLOF',
    'P/LP vs. VUS'
  ),
  mtc = 'fdr',
  alpha = 0.05,
  confint = TRUE,
  xlabel = 'Coefficient Estimate (SD)',
  width = 10,
  height = 8,
  legend_position = 'bottom' #'bottom' for standalone; 'right' when combined with the phenotype distributions
) {
  #MTC is applied across all comparisons x phenotypes (consistent with summary_plotter)
  plot_tb <- input_tb %>%
    mutate(
      adj_p = p.adjust(pval, method = mtc),
      sig = adj_p < alpha,
      lowerCI = estimate - 1.96 * SE,
      upperCI = estimate + 1.96 * SE,
      pheno_short = str_trim(str_match(Label_short, '(^[^()]+)')[, 2]),
      comparison = factor(str_c(group1, ' vs. ', group2), levels = comporder)
    )

  #Symmetric shared x-axis so direction of effect is visually comparable across panels
  x_lim <- max(
    abs(c(plot_tb$lowerCI, plot_tb$upperCI, plot_tb$estimate)),
    na.rm = TRUE
  ) *
    1.05
  y_data_max <- max(-log10(plot_tb$adj_p), -log10(alpha), na.rm = TRUE)
  y_max <- y_data_max * 1.35 #Headroom so the 'Higher in...' text sits clear of the points and their labels
  y_min <- -0.05 * y_max #Small negative buffer so points at adj_p ~ 1 aren't clipped at the bottom

  #Per-panel direction-of-effect annotation (estimate > 0 means higher in group1)
  direction_tb <- plot_tb %>%
    distinct(comparison, group1, group2) %>%
    mutate(
      left_label = str_c('← Higher in ', group2),
      right_label = str_c('Higher in ', group1, ' →')
    )

  #Consistent phenotype colours from colour_hardcoder()
  colour_map <- plot_tb %>%
    distinct(pheno_short, colour) %>%
    drop_na()

  panel_plot <- ggplot(
    plot_tb,
    aes(x = estimate, y = -log10(adj_p), colour = pheno_short)
  ) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey50') +
    geom_hline(yintercept = -log10(alpha), linetype = 'dotted', colour = 'red') +
    {
      if (isTRUE(confint)) {
        geom_linerange(
          aes(xmin = lowerCI, xmax = upperCI),
          alpha = 0.4,
          linewidth = 0.6
        )
      }
    } +
    geom_point(aes(shape = sig), size = 2.8, stroke = 1) +
    ggrepel::geom_text_repel(
      data = filter(plot_tb, sig),
      aes(label = pheno_short),
      size = 3.5,
      show.legend = FALSE,
      min.segment.length = 0,
      box.padding = 0.4,
      ylim = c(NA, y_data_max * 1.15) #Keep repelled labels below the 'Higher in...' text
    ) +
    geom_text(
      data = direction_tb,
      aes(x = -x_lim, y = y_max, label = left_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      colour = 'grey30'
    ) +
    geom_text(
      data = direction_tb,
      aes(x = x_lim, y = y_max, label = right_label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 1,
      size = 3.2,
      colour = 'grey30'
    ) +
    facet_wrap(~comparison, nrow = 2, drop = FALSE) + #Fixed (shared) x and y scales across panels
    scale_colour_manual(
      values = set_names(colour_map$colour, colour_map$pheno_short)
    ) +
    scale_shape_manual(
      values = c(`TRUE` = 16, `FALSE` = 1),
      labels = c(
        `TRUE` = paste0('FDR < ', alpha),
        `FALSE` = paste0('FDR ≥ ', alpha)
      )
    ) +
    scale_x_continuous(limits = c(-x_lim, x_lim), n.breaks = 8) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = function(l) {
        b <- scales::breaks_pretty()(c(0, y_data_max))
        b[b >= 0]
      }, #No ticks/breaks in the negative buffer or the text headroom
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = xlabel,
      y = '-log10(FDR-adjusted p-value)',
      colour = 'Phenotype',
      shape = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_rect(fill = 'grey92'),
      strip.text = element_text(face = 'bold'),
      panel.grid.minor = element_blank(),
      legend.position = legend_position,
      legend.box = 'vertical'
    ) +
    {
      #Two-row legend when at the bottom; single column when at the side
      if (legend_position %in% c('bottom', 'top')) {
        guides(
          colour = guide_legend(nrow = 2, override.aes = list(shape = 16))
        )
      } else {
        guides(
          colour = guide_legend(ncol = 1, override.aes = list(shape = 16))
        )
      }
    }

  ggsave(
    file.path(output_path, file_name),
    panel_plot,
    width = width,
    height = height,
    dpi = 300
  )

  return(panel_plot)
}

#Shared data prep for sig_density_plotter() / sig_violin_plotter(): significant comparison x phenotype pairs,
#their raw values in long format, and per-phenotype facet labels. Returns NULL if nothing is significant
sig_findings_prepper <- function(
  summary_tb,
  comparison_tbs_list,
  param_list,
  comporder,
  mtc,
  alpha
) {
  #Significant comparison x phenotype pairs (MTC across all tests, consistent with summary_plotter)
  sig_tb <- summary_tb %>%
    mutate(adj_p = p.adjust(pval, method = mtc)) %>%
    filter(adj_p < alpha) %>%
    mutate(comparison = str_c(group1, ' vs. ', group2)) #Built from group1/group2 to match the label rebuilt from the comparison_tbs below

  if (nrow(sig_tb) == 0) {
    return(NULL)
  }

  sig_phenos <- unique(sig_tb$pheno)

  #Raw values in long format; comparison label rebuilt from group factor levels (levels = c(group2, group1))
  raw_long <- map_dfr(comparison_tbs_list, function(tb) {
    lv <- levels(tb$group)
    tb %>%
      mutate(
        group = as.character(group),
        comparison = str_c(lv[2], ' vs. ', lv[1])
      ) %>%
      select(HCR_IDs, group, comparison, any_of(sig_phenos)) %>%
      pivot_longer(any_of(sig_phenos), names_to = 'pheno', values_to = 'value')
  }) %>%
    filter(!is.na(value)) %>%
    semi_join(sig_tb, by = c('comparison', 'pheno'))

  #Sub-panel titles: native-unit phenotype label only
  facet_tb <- sig_tb %>%
    mutate(
      facet_label = str_wrap(
        map_chr(pheno, ~ param_list[[.x]][2]),
        35
      )
    ) %>%
    select(comparison, pheno, group1, group2, facet_label)

  list(
    sig_tb = sig_tb,
    raw_long = raw_long,
    facet_tb = facet_tb,
    comps_present = intersect(comporder, unique(sig_tb$comparison))
  )
}

#Per-comparison long tibble (with facet labels) and named group colour vector, shared by the density/violin plotters
sig_comparison_data <- function(prep, comp, palette, default_palette) {
  comp_facets <- filter(prep$facet_tb, comparison == comp)
  g1 <- comp_facets$group1[1]
  g2 <- comp_facets$group2[1]
  comp_palette <- if (!is.null(palette[[comp]])) {
    palette[[comp]]
  } else {
    default_palette
  }

  comp_tb <- prep$raw_long %>%
    filter(comparison == comp) %>%
    left_join(
      select(comp_facets, comparison, pheno, facet_label),
      by = c('comparison', 'pheno')
    ) %>%
    mutate(group = factor(group, levels = c(g1, g2)))

  list(comp_tb = comp_tb, group_cols = set_names(comp_palette, c(g1, g2)))
}

#This plots the raw (non-transformed) phenotype distributions as overlapping histograms across both subgroups
#for every significant comparison x phenotype finding; one row per comparison, one sub-panel per phenotype
sig_density_plotter <- function(
  summary_tb,
  comparison_tbs_list, #List of comparison tibbles each with a `group` factor (levels = c(group2, group1))
  param_list,
  output_path,
  file_name = 'Significant_findings_raw_histogram_panelplot.png',
  comporder = c(
    'Thick vs. Thin',
    'MYBPC3 vs. MYH7',
    'MYBPC3_pLOF vs. MYBPC3_nonLOF',
    'P/LP vs. VUS'
  ),
  mtc = 'fdr',
  alpha = 0.05,
  palette = list(
    #Distinct (group1, group2) colour pair per comparison
    'Thick vs. Thin' = c('#1B9E77', '#D95F02'),
    'MYBPC3 vs. MYH7' = c('#7570B3', '#E6AB02'),
    'MYBPC3_pLOF vs. MYBPC3_nonLOF' = c('#E7298A', '#66A61E'),
    'P/LP vs. VUS' = c('#2166AC', '#B2182B')
  ),
  default_palette = c('grey30', 'grey70'),
  bins = 30,
  hist_scale = c('density', 'count'), #'density' = each group's histogram integrates to 1; 'count' = raw counts
  width = 12,
  height = NULL
) {
  hist_scale <- match.arg(hist_scale)

  prep <- sig_findings_prepper(
    summary_tb,
    comparison_tbs_list,
    param_list,
    comporder,
    mtc,
    alpha
  )

  if (is.null(prep)) {
    message(
      'No significant findings at the chosen threshold; no density plot produced.'
    )
    return(invisible(NULL))
  }

  comp_plots <- map(prep$comps_present, function(comp) {
    comp_data <- sig_comparison_data(prep, comp, palette, default_palette)
    comp_tb <- comp_data$comp_tb
    group_cols <- comp_data$group_cols

    median_tb <- comp_tb %>%
      group_by(facet_label, group) %>%
      summarise(med = median(value), .groups = 'drop')

    ggplot(comp_tb, aes(x = value, fill = group, colour = group)) +
      {
        #Overlapping histograms; density-scaled by default so groups of very different size remain comparable
        if (hist_scale == 'density') {
          geom_histogram(
            aes(y = after_stat(density)),
            position = 'identity',
            bins = bins,
            alpha = 0.4,
            linewidth = 0.3
          )
        } else {
          geom_histogram(
            position = 'identity',
            bins = bins,
            alpha = 0.4,
            linewidth = 0.3
          )
        }
      } +
      geom_vline(
        data = median_tb,
        aes(xintercept = med, colour = group),
        linetype = 'dashed',
        linewidth = 0.7
      ) +
      facet_wrap(~facet_label, scales = 'free', nrow = 1) +
      scale_fill_manual(values = group_cols) +
      scale_colour_manual(values = group_cols) +
      labs(
        title = comp,
        x = 'Raw phenotype value',
        y = ifelse(hist_scale == 'density', 'Density', 'Count'),
        fill = NULL,
        colour = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = 'bold'),
        strip.background = element_blank(), #Plain text phenotype label, no box
        strip.text = element_text(size = 10, face = 'bold'),
        panel.grid.minor = element_blank(),
        legend.position = 'right'
      )
  })

  combined_plot <- patchwork::wrap_plots(comp_plots, ncol = 1) +
    patchwork::plot_annotation(tag_levels = 'A')

  if (is.null(height)) {
    height <- 3.8 * length(comp_plots)
  }

  ggsave(
    file.path(output_path, file_name),
    combined_plot,
    width = width,
    height = height,
    dpi = 300
  )

  return(combined_plot)
}

#This plots the raw (non-transformed) phenotype distributions as side-by-side violins (with inner boxplots) across both subgroups
#for every significant comparison x phenotype finding; one row per comparison, one sub-panel per phenotype
sig_violin_plotter <- function(
  summary_tb,
  comparison_tbs_list, #List of comparison tibbles each with a `group` factor (levels = c(group2, group1))
  param_list,
  output_path,
  file_name = 'Significant_findings_raw_violin_panelplot.png',
  comporder = c(
    'Thick vs. Thin',
    'MYBPC3 vs. MYH7',
    'MYBPC3_pLOF vs. MYBPC3_nonLOF',
    'P/LP vs. VUS'
  ),
  mtc = 'fdr',
  alpha = 0.05,
  palette = list(
    #Distinct (group1, group2) colour pair per comparison
    'Thick vs. Thin' = c('#1B9E77', '#D95F02'),
    'MYBPC3 vs. MYH7' = c('#7570B3', '#E6AB02'),
    'MYBPC3_pLOF vs. MYBPC3_nonLOF' = c('#E7298A', '#66A61E'),
    'P/LP vs. VUS' = c('#2166AC', '#B2182B')
  ),
  default_palette = c('grey30', 'grey70'),
  show_points = FALSE, #Overlay jittered raw points (can be busy for large groups)
  trim = TRUE, #TRUE trims violins to the observed data range
  width = 12,
  height = NULL
) {
  prep <- sig_findings_prepper(
    summary_tb,
    comparison_tbs_list,
    param_list,
    comporder,
    mtc,
    alpha
  )

  if (is.null(prep)) {
    message(
      'No significant findings at the chosen threshold; no violin plot produced.'
    )
    return(invisible(NULL))
  }

  comp_plots <- map(prep$comps_present, function(comp) {
    comp_data <- sig_comparison_data(prep, comp, palette, default_palette)
    comp_tb <- comp_data$comp_tb
    group_cols <- comp_data$group_cols

    ggplot(comp_tb, aes(x = group, y = value, fill = group, colour = group)) +
      geom_violin(
        alpha = 0.4,
        linewidth = 0.4,
        trim = trim,
        scale = 'width' #Equal max width per group so groups of very different size remain comparable
      ) +
      {
        if (isTRUE(show_points)) {
          geom_jitter(width = 0.12, height = 0, size = 0.5, alpha = 0.3)
        }
      } +
      geom_boxplot(
        width = 0.15,
        fill = 'white',
        outlier.shape = NA,
        linewidth = 0.5
      ) + #Box centre line = median
      facet_wrap(~facet_label, scales = 'free_y', nrow = 1) +
      scale_fill_manual(values = group_cols) +
      scale_colour_manual(values = group_cols) +
      labs(
        title = comp,
        x = NULL,
        y = 'Raw phenotype value',
        fill = NULL,
        colour = NULL
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(face = 'bold'),
        strip.background = element_blank(), #Plain text phenotype label, no box
        strip.text = element_text(size = 10, face = 'bold'),
        panel.grid.minor = element_blank(),
        legend.position = 'right'
      )
  })

  combined_plot <- patchwork::wrap_plots(comp_plots, ncol = 1) +
    patchwork::plot_annotation(tag_levels = 'A')

  if (is.null(height)) {
    height <- 3.8 * length(comp_plots)
  }

  ggsave(
    file.path(output_path, file_name),
    combined_plot,
    width = width,
    height = height,
    dpi = 300
  )

  return(combined_plot)
}

#This combines the panelled summary scatter (top) with the raw phenotype distributions for significant findings (underneath)
#into a single publication figure. Panels are tagged sequentially (A = summary scatter; B onwards = one row per comparison)
#dist_type = 'histogram' uses sig_density_plotter() (overlapping histograms); 'violin' uses sig_violin_plotter()
combined_figure_plotter <- function(
  summary_tb,
  comparison_tbs_list,
  param_list,
  output_path,
  dist_type = c('histogram', 'violin'),
  file_name = NULL, #NULL -> 'LinearRegression_summary_panelplot_with_<dist_type>s.png'
  scatter_file_name = 'LinearRegression_summary_panelplot_fdr_legendright.png',
  density_file_name = NULL, #NULL -> 'Significant_findings_raw_<dist_type>_panelplot.png'
  comporder = c(
    'Thick vs. Thin',
    'MYBPC3 vs. MYH7',
    'MYBPC3_pLOF vs. MYBPC3_nonLOF',
    'P/LP vs. VUS'
  ),
  palette = NULL, #NULL uses the distribution plotter's default per-comparison palette
  mtc = 'fdr',
  alpha = 0.05,
  scatter_height = 8, #Height (in) allotted to the summary scatter
  density_row_height = 3.8, #Height (in) allotted to each distribution row
  width = 12,
  ... #Passed to the distribution plotter (e.g. bins/hist_scale for histogram; show_points/trim for violin)
) {
  dist_type <- match.arg(dist_type)
  if (is.null(file_name)) {
    file_name <- str_c(
      'LinearRegression_summary_panelplot_with_',
      dist_type,
      's.png'
    )
  }
  if (is.null(density_file_name)) {
    density_file_name <- str_c(
      'Significant_findings_raw_',
      dist_type,
      '_panelplot.png'
    )
  }

  #Summary scatter with legend on the right (saved under a separate name so the bottom-legend standalone PNG isn't overwritten)
  scatter_plot <- summary_panel_plotter(
    summary_tb,
    output_path = output_path,
    file_name = scatter_file_name,
    comporder = comporder,
    mtc = mtc,
    alpha = alpha,
    legend_position = 'right'
  )

  #Distribution rows for significant findings (re-saves its own standalone PNG as a side effect)
  density_args <- c(
    list(
      summary_tb = summary_tb,
      comparison_tbs_list = comparison_tbs_list,
      param_list = param_list,
      output_path = output_path,
      file_name = density_file_name,
      comporder = comporder,
      mtc = mtc,
      alpha = alpha
    ),
    list(...)
  )
  if (!is.null(palette)) {
    density_args$palette <- palette
  }
  dist_plotter <- switch(
    dist_type,
    histogram = sig_density_plotter,
    violin = sig_violin_plotter
  )
  density_plot <- do.call(dist_plotter, density_args)

  if (is.null(density_plot)) {
    message(
      'No significant findings; combined figure contains the summary scatter only.'
    )
    ggsave(
      file.path(output_path, file_name),
      scatter_plot,
      width = width,
      height = scatter_height,
      dpi = 300
    )
    return(scatter_plot)
  }

  #Number of distribution rows = number of comparisons (in comporder) with >=1 significant finding (same MTC as the plotters)
  n_density_rows <- summary_tb %>%
    mutate(adj_p = p.adjust(pval, method = mtc)) %>%
    filter(adj_p < alpha) %>%
    mutate(comparison = str_c(group1, ' vs. ', group2)) %>%
    filter(comparison %in% comporder) %>%
    distinct(comparison) %>%
    nrow()

  #Nested patchwork tags continue sequentially: A = scatter, B.. = distribution rows
  combined_plot <- patchwork::wrap_plots(
    scatter_plot,
    density_plot,
    ncol = 1,
    heights = c(scatter_height, density_row_height * n_density_rows)
  ) +
    patchwork::plot_annotation(tag_levels = 'A') &
    theme(plot.tag = element_text(face = 'bold', size = 16))

  ggsave(
    file.path(output_path, file_name),
    combined_plot,
    width = width,
    height = scatter_height + density_row_height * n_density_rows,
    dpi = 300,
    limitsize = FALSE
  )

  return(combined_plot)
}
