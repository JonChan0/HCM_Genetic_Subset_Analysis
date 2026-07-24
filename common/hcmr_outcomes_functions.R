data_preprocessor <- function(
  outcomes_df,
  vars_of_interest,
  hcmr_prs = NULL,
  sarcomere_status = c(T, F),
  prs_scale = T,
  hcmr_pcs = NULL,
  log_transform_vars = NULL,
  hcmr_rarevar_class = NULL
) {
  output_tb <- outcomes_df

  if (!is.null(hcmr_prs)) {
    # expects a named list of PRS tibbles
    if (!is.list(hcmr_prs)) {
      stop("hcmr_prs must be a list of tibbles (or NULL).")
    }

    nm <- names(hcmr_prs)
    if (is.null(nm) || any(!nzchar(nm))) {
      stop(
        "hcmr_prs must be a *named* list. Names will be used as new column names."
      )
    }

    for (prs_name in nm) {
      prs_tb <- hcmr_prs[[prs_name]]

      output_tb <- output_tb |>
        left_join(
          dplyr::select(prs_tb, HCR_IDs = IID, !!prs_name := SCORE1_AVG),
          by = "HCR_IDs"
        )

      if (isTRUE(prs_scale)) {
        output_tb <- output_tb |>
          dplyr::mutate(!!prs_name := as.numeric(scale(.data[[prs_name]])[, 1]))
      }
    }
  }

  output_tb <- output_tb %>%
    select(
      HCR_IDs,
      age,
      sex = gender,
      bsa,
      any_of(vars_of_interest),
      starts_with('cv_prs'),
      sarcomere
    )

  if (!is.null(hcmr_rarevar_class)) {
    output_tb <- output_tb %>%
      left_join(hcmr_rarevar_class, by = 'HCR_IDs') %>%
      mutate(
        sarcomere_binned = ifelse(
          is.na(sarcomere_binned),
          'Sarcomere-negative',
          sarcomere_binned
        )
      )
  }

  if (!is.null(hcmr_pcs)) {
    output_tb <- output_tb %>%
      left_join(select(hcmr_pcs, HCR_IDs = IID, PC1:PC10))
  }

  if (!is.null(log_transform_vars)) {
    #Apply log1p transformation to specified variables
    vars_to_log <- intersect(log_transform_vars, names(output_tb))
    if (length(vars_to_log) > 0) {
      output_tb <- output_tb %>%
        mutate(across(all_of(vars_to_log), ~ log1p(.)))
    }
  }

  return(output_tb)
}

# Z-score `predictors` in `df` using the mean/SD fitted on a reference dataset
# (`ref_df`, e.g. the original non-imputed data). This keeps the scaling constant
# across multiply-imputed datasets and identical to the reference's standardisation.
# Mean/SD are computed on observed values (na.rm = TRUE) of the reference.
scale_predictors_fixed <- function(df, predictors, ref_df) {
  preds <- intersect(predictors, intersect(colnames(df), colnames(ref_df)))
  for (p in preds) {
    if (is.numeric(ref_df[[p]]) && is.numeric(df[[p]])) {
      mu <- mean(ref_df[[p]], na.rm = TRUE)
      sdv <- sd(ref_df[[p]], na.rm = TRUE)
      if (is.finite(sdv) && sdv > 0) {
        df[[p]] <- (df[[p]] - mu) / sdv
      }
    }
  }
  df
}

# Helper to get predictor type
get_pred_type <- function(var) {
  if (var %in% continuous_predictors) "continuous" else "categorical"
}


# Cox regression function which takes in the data, either a continuous or categorical predictor variable,
# covariates and specified indicator variable name and time_to_event_or_censoring variable name and
# returns an output tibble of the summary statistics from the analysis for the specified predictor and
# the covariates as well
cox_function <- function(
  input_data,
  predictor_name,
  predictor_type,
  covars,
  indicator_varname,
  time_varname,
  output_dir,
  suffix,
  entry_time_varname = NULL, # name of entry-time variable for left-truncation/delayed entry (e.g. age at recruitment). NULL = right-censored only.
  time_dependent_predictor = NULL, # name(s) of variable(s) (predictor or covariate) to apply tt() to
  time_dependent_type = c("identity", "log", "sqrt"), # transformation for time-varying effect
  age0 = 60, # baseline time at which to report main effect (centering for tt())
  return_model = FALSE,
  run_PHassumption_check = TRUE
) {
  library(survival)

  stopifnot(is.data.frame(input_data))
  stopifnot(is.character(predictor_name), length(predictor_name) == 1)
  stopifnot(is.character(predictor_type), length(predictor_type) == 1)
  stopifnot(is.character(indicator_varname), length(indicator_varname) == 1)
  stopifnot(is.character(time_varname), length(time_varname) == 1)
  stopifnot(is.logical(return_model), length(return_model) == 1)
  stopifnot(
    is.logical(run_PHassumption_check),
    length(run_PHassumption_check) == 1
  )

  predictor_type <- tolower(predictor_type)
  if (!predictor_type %in% c("continuous", "categorical")) {
    stop("predictor_type must be one of: 'continuous', 'categorical'.")
  }

  if (!predictor_name %in% names(input_data)) {
    stop("predictor_name not found in input_data.")
  }
  if (!indicator_varname %in% names(input_data)) {
    stop("indicator_varname not found in input_data.")
  }
  if (!time_varname %in% names(input_data)) {
    stop("time_varname not found in input_data.")
  }
  if (!is.null(entry_time_varname)) {
    if (!is.character(entry_time_varname) || length(entry_time_varname) != 1) {
      stop("entry_time_varname must be NULL or a single column name.")
    }
    if (!entry_time_varname %in% names(input_data)) {
      stop("entry_time_varname not found in input_data.")
    }
    if (identical(entry_time_varname, time_varname)) {
      stop("entry_time_varname must differ from time_varname.")
    }
  }

  if (!is.null(covars)) {
    if (!is.character(covars)) {
      stop("covars must be a character vector of column names (or NULL).")
    }
    missing_covars <- setdiff(covars, names(input_data))
    if (length(missing_covars) > 0) {
      stop(
        "Covariates not found in input_data: ",
        paste(missing_covars, collapse = ", ")
      )
    }
  } else {
    covars <- character(0)
  }

  # ---- Validate time-dependent spec (now supports vector) ----
  time_dependent_type <- match.arg(time_dependent_type)

  if (!is.null(time_dependent_predictor)) {
    if (!is.character(time_dependent_predictor)) {
      stop(
        "time_dependent_predictor must be NULL or a character vector of variable names."
      )
    }
    time_dependent_predictor <- unique(time_dependent_predictor)
    time_dependent_predictor <- time_dependent_predictor[nzchar(
      time_dependent_predictor
    )]
    if (length(time_dependent_predictor) == 0) time_dependent_predictor <- NULL
  }

  if (!is.null(time_dependent_predictor)) {
    missing_td <- setdiff(time_dependent_predictor, names(input_data))
    if (length(missing_td) > 0) {
      stop(
        "time_dependent_predictor not found in input_data: ",
        paste(missing_td, collapse = ", ")
      )
    }

    allowed <- unique(c(predictor_name, covars))
    not_allowed <- setdiff(time_dependent_predictor, allowed)
    if (length(not_allowed) > 0) {
      stop(
        "time_dependent_predictor must be either predictor_name or one of covars. Offenders: ",
        paste(not_allowed, collapse = ", ")
      )
    }

    # ensure time-dependent vars appear as main effects in RHS
    to_add <- setdiff(time_dependent_predictor, c(predictor_name, covars))
    if (length(to_add) > 0) covars <- unique(c(covars, to_add))
  }

  # Keep only needed columns and drop missing values for model variables
  model_vars <- unique(c(
    entry_time_varname,
    time_varname,
    indicator_varname,
    predictor_name,
    covars,
    time_dependent_predictor
  ))
  dat <- input_data |>
    dplyr::select(dplyr::all_of(model_vars)) |>
    dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything())))

  # Left-truncation: drop rows with non-positive risk interval (entry >= exit);
  # these contribute nothing to any risk set and would otherwise error in coxph.
  if (!is.null(entry_time_varname)) {
    bad <- dat[[entry_time_varname]] >= dat[[time_varname]]
    if (any(bad, na.rm = TRUE)) {
      warning(
        sum(bad, na.rm = TRUE),
        " row(s) dropped where ",
        entry_time_varname,
        " >= ",
        time_varname,
        " (non-positive risk interval) when testing ",
        suffix
      )
      dat <- dat[!bad, , drop = FALSE]
    }
  }

  # Basic type handling for main predictor only (as before)
  if (predictor_type == "categorical") {
    dat <- dat |>
      dplyr::mutate(
        !!predictor_name := as.factor(.data[[predictor_name]])
      )
  } else {
    dat <- dat |>
      dplyr::mutate(
        !!predictor_name := as.numeric(.data[[predictor_name]])
      )
  }

  # Handle time-dependent variables: must be numeric; if binary categorical, recode to 0/1
  # (tt() supports numeric x; multi-level categorical not supported here)
  if (!is.null(time_dependent_predictor)) {
    for (v in time_dependent_predictor) {
      x <- dat[[v]]
      is_binary <- FALSE

      if (is.factor(x)) {
        is_binary <- nlevels(x) == 2
        if (is_binary) dat[[v]] <- as.numeric(x == levels(x)[2])
      } else if (is.character(x)) {
        ux <- unique(x[!is.na(x)])
        is_binary <- length(ux) == 2
        if (is_binary) dat[[v]] <- as.numeric(x == ux[2])
      } else if (is.logical(x)) {
        is_binary <- TRUE
        dat[[v]] <- as.numeric(x)
      } else if (is.numeric(x) || is.integer(x)) {
        ux <- sort(unique(x[is.finite(x)]))
        is_binary <- length(ux) == 2 && all(ux %in% c(0, 1))
        dat[[v]] <- as.numeric(x)
      }

      if (!(is.numeric(dat[[v]]) || is.integer(dat[[v]]))) {
        stop(
          "time_dependent_predictor '",
          v,
          "' must be numeric or binary (2-level factor/character or TRUE/FALSE or 0/1)."
        )
      }
    }
  }

  # tt() function
  tt_fn <- NULL
  if (!is.null(time_dependent_predictor)) {
    tt_fn <- switch(
      time_dependent_type,
      log = function(x, t, ...) x * (log(t) - log(age0)),
      sqrt = function(x, t, ...) x * (sqrt(t) - sqrt(age0)),
      identity = function(x, t, ...) x * (t - age0)
    )
  }

  # Build RHS with tt(var) for each nominated variable
  rhs_terms <- c(predictor_name, covars)
  rhs_terms <- rhs_terms[rhs_terms != ""]
  rhs <- if (length(rhs_terms) == 0) "1" else paste(rhs_terms, collapse = " + ")

  if (!is.null(time_dependent_predictor)) {
    rhs <- paste0(
      rhs,
      " + ",
      paste0("tt(", time_dependent_predictor, ")", collapse = " + ")
    )
  }

  surv_lhs <- if (is.null(entry_time_varname)) {
    paste0("Surv(", time_varname, ", ", indicator_varname, ")")
  } else {
    # Delayed-entry / left-truncated: Surv(entry, exit, event)
    paste0(
      "Surv(",
      entry_time_varname,
      ", ",
      time_varname,
      ", ",
      indicator_varname,
      ")"
    )
  }

  fml <- stats::as.formula(paste0(surv_lhs, " ~ ", rhs))

  fit <- if (is.null(time_dependent_predictor)) {
    survival::coxph(fml, data = dat)
  } else {
    survival::coxph(fml, data = dat, tt = tt_fn)
  }

  # ---- Optional PH assumption check (cox.zph + BH correction) ----
  zph_tbl <- NULL
  ph_failed_vars <- character(0)

  if (isTRUE(run_PHassumption_check)) {
    # cox.zph is not defined for models that already include tt() time-varying effects
    if (is.null(time_dependent_predictor)) {
      zph <- survival::cox.zph(fit)

      zph_tbl <- as.data.frame(zph$table) |>
        tibble::rownames_to_column("term") |>
        dplyr::rename(
          chisq = chisq,
          df = df,
          p.value = p
        ) |>
        dplyr::filter(term != "GLOBAL") |>
        dplyr::mutate(
          adj_pvalue_zph = stats::p.adjust(p.value, method = "holm")
        )

      flagged <- zph_tbl |>
        dplyr::filter(is.finite(adj_pvalue_zph), adj_pvalue_zph < 0.05)

      if (nrow(flagged) > 0) {
        ph_failed_vars <- flagged$term
        apply(flagged, 1, function(r) {
          print(paste0(
            "Statistically significant deviation from PH assumption detected for ",
            r[["term"]],
            " with adjusted p-value ",
            r[["adj_pvalue_zph"]],
            " when testing ",
            suffix
          ))
        })
      }
    }
  }

  if (isTRUE(return_model)) {
    attr(fit, "ph_failed_vars") <- ph_failed_vars
    return(fit)
  }

  tidy_tbl <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) |>
    dplyr::transmute(
      term = term,
      estimate = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value
    ) |>
    dplyr::rename(
      hr = estimate,
      hr_low = conf.low,
      hr_high = conf.high
    )

  # Add meta + sample size info
  out <- tidy_tbl |>
    dplyr::mutate(
      predictor = predictor_name,
      predictor_type = predictor_type,
      time_dependent_predictor = if (is.null(time_dependent_predictor)) {
        NA_character_
      } else {
        paste(time_dependent_predictor, collapse = ",")
      },
      n = stats::nobs(fit),
      n_events = sum(dat[[indicator_varname]] == 1, na.rm = TRUE),
      outcome_entry_var = if (is.null(entry_time_varname)) {
        NA_character_
      } else {
        entry_time_varname
      },
      outcome_time_var = time_varname,
      outcome_event_var = indicator_varname,
      .before = 1
    ) |>
    dplyr::mutate(composite = stringr::str_remove(suffix, "^_")) # remove leading underscore for composite column

  attr(out, "ph_failed_vars") <- ph_failed_vars

  out
}

# Wrapper to run Cox and handle PH failure automatically via time-dependent modeling
run_adaptive_cox <- function(
  data,
  predictor,
  predictor_type,
  covariates,
  outcome,
  time,
  suffix,
  output_dir,
  entry_time = NULL, # entry-time variable for left-truncation/delayed entry; passed through to cox_function
  run_PHassumption_check = TRUE
) {
  # Initial run with PH check
  res_tibble <- cox_function(
    input_data = data,
    predictor_name = predictor,
    predictor_type = predictor_type,
    covars = covariates,
    indicator_varname = outcome,
    time_varname = time,
    entry_time_varname = entry_time,
    output_dir = output_dir,
    suffix = suffix,
    run_PHassumption_check = run_PHassumption_check
  )

  # Check for PH failure
  if (run_PHassumption_check) {
    ph_failed_vars <- attr(res_tibble, "ph_failed_vars")

    if (length(ph_failed_vars) > 0) {
      message(paste(
        "PH assumption failed for:",
        paste(ph_failed_vars, collapse = ", "),
        "in model:",
        predictor,
        suffix,
        "- Rerunning with time-dependent covariates."
      ))

      # Determine age to center by. Gracefully handle if 'age' doesn't exist.
      age_center <- if ("age" %in% names(data)) {
        median(data[['age']], na.rm = TRUE)
      } else {
        0
      }

      res_tibble <- cox_function(
        input_data = data,
        predictor_name = predictor,
        predictor_type = predictor_type,
        covars = covariates,
        indicator_varname = outcome,
        time_varname = time,
        entry_time_varname = entry_time,
        output_dir = output_dir,
        suffix = paste0(suffix, "_tt"),
        time_dependent_predictor = ph_failed_vars,
        age0 = age_center,
        time_dependent_type = "identity",
        run_PHassumption_check = FALSE
      )
    }
  }

  return(res_tibble)
}

# Diagnostics for Cox regression models:
# - Proportional hazards: cox.zph test + scaled Schoenfeld residual plots
# - Functional form (linearity): martingale residual plots for continuous terms
# - Influence: deviance residuals, dfbeta
# Saves plots + a small list of test tibbles to `output_dir`

cox_diagnostics <- function(
  fit,
  output_dir,
  file_stub = "cox_diag",
  suffix = "",
  ph_transform = "km",
  ph_terms = "all",
  martingale_smoother = TRUE,
  return_objects = TRUE
) {
  stopifnot(inherits(fit, "coxph"))
  stopifnot(is.character(output_dir), length(output_dir) == 1)
  stopifnot(is.character(file_stub), length(file_stub) == 1)
  stopifnot(is.character(suffix), length(suffix) == 1)
  stopifnot(is.character(ph_transform), length(ph_transform) == 1)
  stopifnot(is.logical(martingale_smoother), length(martingale_smoother) == 1)
  stopifnot(is.logical(return_objects), length(return_objects) == 1)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  safe_stub <- gsub("[^A-Za-z0-9_\\-]+", "_", paste0(file_stub, suffix))
  path <- function(name) file.path(output_dir, paste0(safe_stub, "_", name))

  # ---- 1) Proportional hazards assumption (Schoenfeld residuals) ----
  zph <- survival::cox.zph(fit, transform = ph_transform)

  # Save test table
  zph_tbl <- as.data.frame(zph$table) |>
    tibble::rownames_to_column("term") |>
    dplyr::rename(
      chisq = chisq,
      df = df,
      p.value = p
    )

  readr::write_csv(zph_tbl, paste0(path("ph_test"), ".csv"))

  # Save Schoenfeld plots (base graphics from survival)
  # If ph_terms = "all", plot all covariates; otherwise allow a character vector of term names.
  ph_terms_to_plot <- ph_terms
  if (identical(ph_terms_to_plot, "all")) {
    ph_terms_to_plot <- seq_len(nrow(zph$table) - 1) # exclude global row (last)
  }

  fn_ph <- paste0(path("ph_schoenfeld"), ".png")
  ragg::agg_png(fn_ph, width = 12, height = 4, units = "in", res = 600)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(
    {
      try(graphics::par(op), silent = TRUE)
      try(grDevices::dev.off(), silent = TRUE)
    },
    add = TRUE
  )

  # Multi-panel layout, one per term (if many, wrap in columns)
  n_panels <- if (identical(ph_terms, "all")) {
    nrow(zph$table) - 1
  } else {
    length(ph_terms_to_plot)
  }
  ncol <- min(3, max(1, ceiling(sqrt(n_panels))))
  nrow <- ceiling(n_panels / ncol)
  graphics::par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))

  plot(zph, resid = TRUE, se = TRUE, var = ph_terms_to_plot)

  # Close device explicitly here
  grDevices::dev.off()

  # ---- 2) Linearity / functional form (Martingale residuals) ----
  # Guidance: for continuous covariates, martingale residuals vs covariate should be ~ linear around 0.
  mf <- stats::model.frame(fit)
  term_labels <- attr(stats::terms(fit), "term.labels")

  # identify "continuous-looking" covariates available in model frame
  cont_terms <- term_labels[
    term_labels %in% names(mf) & vapply(mf[term_labels], is.numeric, logical(1))
  ]

  mart <- stats::residuals(fit, type = "martingale")
  mart_tbl <- tibble::tibble(.martingale = as.numeric(mart))

  if (length(cont_terms) > 0) {
    fn_mart <- paste0(path("linearity_martingale"), ".png")
    ragg::agg_png(fn_mart, width = 12, height = 4, units = "in", res = 600)

    n_panels <- length(cont_terms)
    ncol <- min(3, max(1, ceiling(sqrt(n_panels))))
    nrow <- ceiling(n_panels / ncol)
    graphics::par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))

    for (v in cont_terms) {
      x <- mf[[v]]
      ok <- is.finite(x) & is.finite(mart)
      graphics::plot(
        x[ok],
        mart[ok],
        xlab = v,
        ylab = "Martingale residual",
        main = paste0("Linearity check: ", v),
        pch = 16,
        col = grDevices::rgb(0, 0, 0, 0.35)
      )
      graphics::abline(h = 0, lty = 2, col = "grey50")

      if (isTRUE(martingale_smoother)) {
        graphics::lines(
          stats::lowess(x[ok], mart[ok], f = 2 / 3),
          col = "firebrick",
          lwd = 2
        )
      }
    }

    grDevices::dev.off()
  }

  # ---- 3) Basic influence / residual diagnostics ----
  dev <- stats::residuals(fit, type = "deviance")
  fn_dev <- paste0(path("residuals_deviance"), ".png")
  ragg::agg_png(fn_dev, width = 8, height = 6, units = "in", res = 600)
  graphics::par(mar = c(4, 4, 2, 1))
  graphics::plot(
    dev,
    ylab = "Deviance residual",
    xlab = "Observation index",
    main = "Deviance residuals (index plot)",
    pch = 16,
    col = grDevices::rgb(0, 0, 0, 0.35)
  )
  graphics::abline(h = 0, lty = 2, col = "grey50")
  grDevices::dev.off()

  # DFBETA: one panel per coefficient
  dfb <- stats::residuals(fit, type = "dfbeta")
  if (!is.null(dfb) && ncol(as.matrix(dfb)) > 0) {
    dfb_mat <- as.matrix(dfb)
    fn_dfb <- paste0(path("influence_dfbeta"), ".png")
    ragg::agg_png(fn_dfb, width = 12, height = 4, units = "in", res = 600)

    n_panels <- ncol(dfb_mat)
    ncolp <- min(3, max(1, ceiling(sqrt(n_panels))))
    nrowp <- ceiling(n_panels / ncolp)
    graphics::par(mfrow = c(nrowp, ncolp), mar = c(4, 4, 2, 1))

    for (j in seq_len(ncol(dfb_mat))) {
      nm <- colnames(dfb_mat)[j]
      graphics::plot(
        dfb_mat[, j],
        ylab = "DFBETA",
        xlab = "Observation index",
        main = paste0("DFBETA: ", nm),
        pch = 16,
        col = grDevices::rgb(0, 0, 0, 0.35)
      )
      graphics::abline(h = 0, lty = 2, col = "grey50")
    }

    grDevices::dev.off()
  }

  if (!isTRUE(return_objects)) {
    return(invisible(NULL))
  }

  return(zph_tbl)
}


#Forest plot for Cox regression results
forest_plot_cox <- function(
  input_tbl,
  output_dir,
  fill_by = "Group",
  shape_point = 23,
  group_palette = 'Set2',
  hr_ref = 1,
  x_limits = NULL,
  facet_variable = "composite",
  y_legend_spacing = 10,
  suffix = "",
  ylabel = '',
  legend_position = 'right',
  star_xnudge = 0,
  star_ynudge = 0.1,
  forest_width = 9,
  forest_height = 4
) {
  stopifnot(is.data.frame(input_tbl))

  required_cols <- c(
    "term",
    "n_events",
    "composite",
    "hr",
    "hr_low",
    "hr_high",
    "p.value",
    "Group"
  )
  missing_cols <- setdiff(required_cols, names(input_tbl))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!fill_by %in% names(input_tbl)) {
    stop("fill_by column not found in input_tbl: ", fill_by)
  }
  if (!facet_variable %in% names(input_tbl)) {
    stop("facet_variable column not found in input_tbl: ", facet_variable)
  }

  if ('adj_pvalue' %in% names(input_tbl)) {
    input_tbl <- dplyr::rename(input_tbl, 'adj_p' = 'adj_pvalue')
  }

  has_adj_p <- "adj_p" %in% names(input_tbl)
  p_col <- if (has_adj_p) "adj_p" else "p.value"
  has_custom_annot <- "custom_annot" %in% names(input_tbl)

  dat <- input_tbl |>
    dplyr::mutate(
      term = as.character(term), # Ensure term is character for alphabetical sorting
      composite = as.character(composite),
      Group = as.factor(Group),
      n_events = as.numeric(n_events),
      hr = as.numeric(hr),
      hr_low = as.numeric(hr_low),
      hr_high = as.numeric(hr_high),
      p.value = as.numeric(p.value),
      adj_p = if (has_adj_p) as.numeric(adj_p) else NA_real_,
      .p_for_stars = as.numeric(.data[[p_col]]),
      .stars = if (has_custom_annot) {
        as.character(custom_annot)
      } else {
        dplyr::case_when(
          is.na(.p_for_stars) ~ "",
          .p_for_stars < 0.001 ~ "***",
          .p_for_stars < 0.01 ~ "**",
          .p_for_stars < 0.05 ~ "*",
          TRUE ~ ""
        )
      },
      .fill = .data[[fill_by]],
      .facet = as.character(.data[[facet_variable]])
    ) |>
    dplyr::filter(
      is.finite(hr),
      is.finite(hr_low),
      is.finite(hr_high)
    ) |>
    dplyr::arrange(desc(Group), desc(term)) |> # Order by Group level first, then Alphabetically
    dplyr::mutate(term = fct_inorder(term))

  if (is.null(x_limits)) {
    x_min <- min(dat$hr_low, hr_ref, na.rm = TRUE)
    x_max <- max(dat$hr_high, hr_ref, na.rm = TRUE)
    x_limits <- c(x_min, x_max)
  }

  p <- ggplot2::ggplot(
    dat,
    ggplot2::aes(
      x = hr,
      y = term,
      fill = .fill
    )
  ) +
    ggplot2::geom_vline(
      xintercept = hr_ref,
      linetype = "dashed",
      color = "grey50"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = hr_low, xmax = hr_high),
      height = 0.2,
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      shape = shape_point,
      size = 3,
      color = "black",
      stroke = 0.4
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .stars),
      nudge_x = star_xnudge,
      nudge_y = star_ynudge,
      size = 5,
      vjust = 0,
      show.legend = FALSE
    ) +
    ggplot2::scale_x_log10(limits = x_limits) +
    ggplot2::facet_wrap(ggplot2::vars(.facet), nrow = 1, scales = "fixed") +
    ggplot2::scale_fill_brewer(palette = group_palette) +
    ggplot2::labs(
      x = "Hazard Ratio (log scale)",
      y = ylabel,
      fill = fill_by
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey95"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.key.spacing.y = grid::unit(y_legend_spacing, "points"),
      legend.position = legend_position
    )

  if (!missing(output_dir) && !is.null(output_dir) && nzchar(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    ggplot2::ggsave(
      filename = file.path(
        output_dir,
        str_glue("forest_plot_cox_comp{suffix}.png")
      ),
      plot = p,
      width = forest_width,
      height = forest_height,
      dpi = 600
    )
  }

  p
}
# Kaplan Meier survival curve plotter function
km_survival_curv_plotter <- function(
  input_tb,
  categorical_var,
  status,
  outcome_var,
  plot_output_folder,
  entry_time_varname = NULL,
  ylimits = NULL,
  plot_types = c('pct'),
  conf_int_style = 'ribbon', #Either 'ribbon' or 'step'
  output_suffix = '',
  xlab = 'Age',
  km_width = 12,
  km_height = 6
) {
  # --- Check for necessary packages ---
  if (!requireNamespace("survminer", quietly = TRUE)) {
    stop(
      "Package 'survminer' is needed for this function. Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      "Package 'survival' is needed for this function. Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is needed for this function. Please install it.",
      call. = FALSE
    )
  }

  library(survival)

  # categorical_var, status, and outcome_var are now expected as strings
  print(paste('Plotting Kaplan Meier survival curves for', categorical_var))

  # Check if ylimits is valid
  valid_ylimits <- !is.null(ylimits) &&
    is.list(ylimits) &&
    length(ylimits) == 3 &&
    all(sapply(ylimits, function(y) is.numeric(y) && length(y) == 2))

  if (!is.null(ylimits) && !valid_ylimits) {
    print(
      "Warning: 'ylimits' argument is not NULL but is not a valid list of three 2-element vectors. Ignoring custom ylimits."
    )
  }

  # Construct the survival formula as a formula object.
  # When entry_time_varname is supplied, build a left-truncated (delayed-entry)
  # Surv object: Surv(entry, stop, status). This makes the risk set at each time
  # only include subjects already under observation, which is required for valid
  # age-as-timescale log-rank tests and KM curves.
  if (!is.null(entry_time_varname)) {
    if (is.null(input_tb[[entry_time_varname]])) {
      stop(
        paste0(
          "entry_time_varname '",
          entry_time_varname,
          "' not found in input data."
        ),
        call. = FALSE
      )
    }
    surv_formula <- as.formula(paste0(
      "Surv(",
      entry_time_varname,
      ", ",
      outcome_var,
      ", ",
      status,
      ") ~ ",
      categorical_var
    ))
    print(paste0(
      "Applying left truncation: entry time = '",
      entry_time_varname,
      "', stop time = '",
      outcome_var,
      "'."
    ))
  } else {
    surv_formula <- as.formula(paste0(
      "Surv(",
      outcome_var,
      ", ",
      status,
      ") ~ ",
      categorical_var
    ))
  }

  # Log-rank test via the score (logrank) test of a Cox model. survdiff() does
  # not support counting-process (left-truncated) Surv objects, whereas coxph()
  # does; its score test is the log-rank test and is identical to survdiff() in
  # the non-truncated case. This keeps a single code path for both timescales.
  logrank_pvalue <- function(data) {
    fit <- survival::coxph(surv_formula, data = data)
    unname(summary(fit)$sctest["pvalue"])
  }

  # --- Determine number of groups for pairwise comparison ---
  if (is.null(input_tb[[categorical_var]])) {
    stop(
      paste(
        "Categorical variable '",
        categorical_var,
        "' not found in input data."
      ),
      call. = FALSE
    )
  }
  if (!is.factor(input_tb[[categorical_var]])) {
    input_tb[[categorical_var]] <- as.factor(input_tb[[categorical_var]])
  }
  num_groups <- length(levels(input_tb[[categorical_var]]))

  # Set flag to run manual pairwise tests
  run_manual_pairwise <- num_groups > 2

  if (run_manual_pairwise) {
    print(paste(
      "Found",
      num_groups,
      "groups. Will perform manual pairwise log-rank tests."
    ))

    input_tb[[categorical_var]] <- as_factor(input_tb[[categorical_var]]) #Make categorical variable

    # --- BEGIN NEW FIX for '0,1,2' BUG ---
    # Check if factor levels are "numeric-like" (e.g., "0", "1", "2")
    # This confuses ggplot2's aesthetic mapping (treating discrete as continuous)
    current_levels <- levels(input_tb[[categorical_var]])

    # Check if all non-NA levels can be converted to numbers
    numeric_levels <- suppressWarnings(as.numeric(current_levels))

    if (!any(is.na(numeric_levels))) {
      # All levels are numeric-like. Prepend a string to force discrete treatment.
      print(paste(
        "Numeric-like factor levels detected:",
        paste(current_levels, collapse = ", ")
      ))
      print(
        "Prepending 'Dosage-' to levels to ensure discrete mapping in ggplot2."
      )

      # Use the levels() function to modify the factor levels directly
      levels(input_tb[[categorical_var]]) <- paste0("Dosage-", current_levels)
    }
    # --- END NEW FIX ---

    print(levels(input_tb[[categorical_var]]))
  } else {
    print(paste(
      "Found",
      num_groups,
      "group(s). Will show single overall p-value."
    ))
  }

  # Use an index to iterate through plot_types and ylimits simultaneously
  for (i in 1:length(plot_types)) {
    type <- plot_types[i]

    # Fit the survival model using the constructed formula
    surv_fit <- survival::survfit(surv_formula, data = input_tb)

    # Update the call in the fitted object
    surv_fit$call$data <- quote(input_tb)
    surv_fit$call$formula <- surv_formula

    # Generate the survival plot using ggsurvplot
    plot_obj <- survminer::ggsurvplot(
      fit = surv_fit,
      data = input_tb,
      fun = type,
      pval = F,
      conf.int = T,
      conf.int.style = conf_int_style,
      risk.table = 'abs_pct',
      risk.table.col = 'strata',
      linetype = ifelse(run_manual_pairwise, 1, 'strata'),
      break.time.by = ifelse(str_detect(xlab, 'Age'), 10, 2),
      xlab = xlab,
      ncensor.plot = F,
      font.x = 16,
      font.y = 16,
      font.tickslab = 14,
      font.legend = 14,
      ggtheme = ggplot2::theme_classic()
    )

    # --- Apply custom y-limits if provided ---
    if (valid_ylimits) {
      current_ylim <- ylimits[[i]]
      if (!is.null(current_ylim)) {
        plot_obj$plot <- plot_obj$plot + ggplot2::ylim(current_ylim)
        print(paste(
          "Applied custom y-limits for '",
          type,
          "': c(",
          current_ylim[1],
          ", ",
          current_ylim[2],
          ")",
          sep = ""
        ))
      }
    }

    # --- Truncation-aware risk table ---
    # survminer's risk table ignores delayed entry: it counts every subject as
    # present from time 0, so for an age timescale it reports a full risk set at
    # ages below the minimum entry age. Recompute n at risk as the number under
    # observation at each break time, i.e. #{entry < t <= stop} within each
    # group, so the table is consistent with the left-truncated curve and test.
    # 'abs_pct' renders the 'llabels' column ("n.risk (pct.risk)"). survminer
    # attaches the table data to the geom_text *layer*, not just the plot-level
    # default, so we must rewrite the layer data for the numbers to actually
    # change on the rendered plot.
    if (!is.null(entry_time_varname) && !is.null(plot_obj$table)) {
      entry_times <- input_tb[[entry_time_varname]]
      stop_times <- input_tb[[outcome_var]]
      grp_vec <- as.character(input_tb[[categorical_var]])

      recompute_risk_table <- function(tbl_data) {
        if (
          !is.data.frame(tbl_data) ||
            !all(c("strata", "time", "llabels") %in% names(tbl_data))
        ) {
          return(tbl_data)
        }

        # survfit names strata "<var>=<level>"; strip the prefix to recover the
        # raw group level (no-op if survminer already cleaned it to a bare level).
        strata_chr <- sub(
          paste0("^", categorical_var, "="),
          "",
          as.character(tbl_data$strata)
        )

        tbl_data$n.risk <- mapply(
          function(s, t) {
            in_grp <- grp_vec == s & !is.na(entry_times) & !is.na(stop_times)
            sum(entry_times[in_grp] < t & stop_times[in_grp] >= t)
          },
          strata_chr,
          tbl_data$time
        )

        if ("strata_size" %in% names(tbl_data)) {
          tbl_data$pct.risk <- round(
            tbl_data$n.risk * 100 / tbl_data$strata_size
          )
          tbl_data$llabels <- paste0(
            tbl_data$n.risk,
            " (",
            tbl_data$pct.risk,
            ")"
          )
        } else {
          tbl_data$llabels <- as.character(tbl_data$n.risk)
        }
        tbl_data
      }

      # Plot-level default data and every geom layer that holds its own copy.
      plot_obj$table$data <- recompute_risk_table(plot_obj$table$data)
      for (li in seq_along(plot_obj$table$layers)) {
        plot_obj$table$layers[[li]]$data <- recompute_risk_table(
          plot_obj$table$layers[[li]]$data
        )
      }
      print(
        "Recomputed risk table accounting for left truncation (delayed entry)."
      )
    }

    # --- BEGIN MANUAL PAIRWISE P-VALUE ADDITION ---
    if (run_manual_pairwise) {
      # 1. Get all unique pairwise combinations of group levels
      group_levels <- levels(input_tb[[categorical_var]])
      combinations <- utils::combn(group_levels, 2, simplify = FALSE)

      pval_strings <- c() # To store formatted p-value strings

      # 2. Loop through each pair, subset data, and run log-rank test
      for (pair in combinations) {
        # Subset the data to only include the two groups in the current pair
        subset_tb <- input_tb[input_tb[[categorical_var]] %in% pair, ]

        # IMPORTANT: Refactor the subsetted data to drop unused levels
        # This is critical for the test to run on only two groups
        subset_tb[[categorical_var]] <- factor(subset_tb[[categorical_var]])

        # Perform the log-rank test (Cox score test) on the pair
        p_value <- logrank_pvalue(subset_tb)

        # Format the p-value and group names into a string
        p_val_formatted <- format.pval(p_value, digits = 3, eps = 0.001)
        current_string <- paste(
          pair[1],
          "vs",
          pair[2],
          ": p =",
          p_val_formatted
        )

        # Add the string to our vector
        pval_strings <- c(pval_strings, current_string)
      }

      # 3. Combine all pairwise strings into a single block of text
      final_label_string <- paste(pval_strings, collapse = "\n")

      # 4. Define coordinates for the text block
      x_pos <- max(input_tb[[outcome_var]], na.rm = TRUE) * 0.15

      y_lower_bound <- 0
      y_upper_bound <- 1

      # Determine default upper bound based on plot type
      if (type == 'cumhaz') {
        max_haz <- max(surv_fit$cumhaz, na.rm = TRUE)
        if (is.finite(max_haz) && max_haz > 0) {
          y_upper_bound <- max_haz
        }
        # if max_haz is 0 or Inf, y_upper_bound remains 1 as a fallback
      }

      # Check if valid custom limits for THIS plot type are provided
      if (valid_ylimits && !is.null(ylimits[[i]])) {
        y_lower_bound <- ylimits[[i]][1]
        y_upper_bound <- ylimits[[i]][2]
        print(paste(
          "Using custom y-range for annotation:",
          y_lower_bound,
          y_upper_bound
        ))
      }

      # Place text at 25% up from the bottom of the visible y-range
      y_pos <- y_lower_bound + (y_upper_bound - y_lower_bound) * 0.25

      # 5. Add the annotation text to the plot
      plot_obj$plot <- plot_obj$plot +
        ggplot2::annotate(
          geom = "text",
          x = x_pos,
          y = y_pos,
          label = final_label_string,
          hjust = 0, # Left-align the text block
          vjust = 0, # Bottom-align the text block (anchor is bottom-left)
          size = 4 # Small text size to avoid clutter
        )
      print("Manually added pairwise p-values to plot.")
    } else {
      # compute the overall log-rank p-value (Cox score test)
      p_value <- logrank_pvalue(input_tb)

      # 4. Define coordinates for the text block
      x_pos <- max(input_tb[[outcome_var]], na.rm = TRUE) * 0.15

      y_lower_bound <- 0
      y_upper_bound <- 1

      # Determine default upper bound based on plot type
      if (type == 'cumhaz') {
        max_haz <- max(surv_fit$cumhaz, na.rm = TRUE)
        if (is.finite(max_haz) && max_haz > 0) {
          y_upper_bound <- max_haz
        }
        # if max_haz is 0 or Inf, y_upper_bound remains 1 as a fallback
      }

      # Check if valid custom limits for THIS plot type are provided
      if (valid_ylimits && !is.null(ylimits[[i]])) {
        y_lower_bound <- ylimits[[i]][1]
        y_upper_bound <- ylimits[[i]][2]
        print(paste(
          "Using custom y-range for annotation:",
          y_lower_bound,
          y_upper_bound
        ))
      }

      # Place text at 25% up from the bottom of the visible y-range
      y_pos <- y_lower_bound + (y_upper_bound - y_lower_bound) * 0.25

      # 5. Add the annotation text to the plot
      plot_obj$plot <- plot_obj$plot +
        ggplot2::annotate(
          geom = "text",
          x = x_pos,
          y = y_pos,
          label = paste0("p = ", signif(p_value, 3)),
          hjust = 0, # Left-align the text block
          vjust = 0, # Bottom-align the text block (anchor is bottom-left)
          size = 6 # Small text size to avoid clutter
        )
      print("Added overall p-value to plot.")
    }
    # --- END MANUAL ADDITION ---

    # Create output folder if it does not exist
    if (!dir.exists(plot_output_folder)) {
      dir.create(plot_output_folder, recursive = TRUE)
      message("Created output folder: ", plot_output_folder)
    }

    png_filename <- file.path(
      plot_output_folder,
      paste0(categorical_var, "_kaplan_", type, output_suffix, ".png")
    )

    tryCatch(
      {
        ragg::agg_png(
          filename = png_filename,
          width = km_width,
          height = km_height,
          units = "in",
          res = 600
        )
        print(plot_obj, newpage = FALSE)
        grDevices::dev.off()
        message("Saved PNG: ", png_filename)
      },
      error = function(e) {
        # Clean up device if open
        try(grDevices::dev.off(), silent = TRUE)

        is_aesthetic_error <- grepl(
          "geom to grob|Aesthetics can not vary",
          conditionMessage(e)
        )
        if (is_aesthetic_error) {
          warning(
            "PNG saving failed (aesthetic conflict). Falling back to PDF."
          )
          warning("Original error: ", conditionMessage(e))

          pdf_filename <- file.path(
            plot_output_folder,
            paste0(categorical_var, "_kaplan_", type, output_suffix, ".pdf")
          )

          grDevices::pdf(pdf_filename, width = 12, height = 9)
          print(plot_obj, newpage = FALSE)
          grDevices::dev.off()
          message("Saved PDF (fallback): ", pdf_filename)

          # Remove any broken/empty PNG that may exist
          try(unlink(png_filename), silent = TRUE)
        } else {
          stop(
            "PNG saving failed with an unexpected error: ",
            conditionMessage(e)
          )
        }
      }
    )
  }
  print(paste('Finished plotting for', categorical_var))
}

#This checks the assumptions of Cox regression and outputs diagnostic plots and tables to a specified folder.
check_cox_assumptions <- function(fit, output_folder, data) {
  # --- 2. Check/Create Output Folder ---
  if (!dir.exists(output_folder)) {
    tryCatch(
      {
        dir.create(output_folder, recursive = TRUE)
        message("Created output folder: ", output_folder)
      },
      warning = function(w) {
        stop(
          "Could not create output folder: ",
          output_folder,
          "\nWarning: ",
          w$message,
          call. = FALSE
        )
      },
      error = function(e) {
        stop(
          "Could not create output folder: ",
          output_folder,
          "\nError: ",
          e$message,
          call. = FALSE
        )
      }
    )
  }

  # --- 3. Proportional Hazards (PH) Assumption ---
  message("--- Checking Proportional Hazards (PH) Assumption ---")
  message("Testing using Schoenfeld residuals (cox.zph)...")

  # Calculate Schoenfeld residuals
  tryCatch(
    {
      zph_result <- survival::cox.zph(fit)

      # Print statistical test results
      message(
        "Test statistics and p-values (a p-value < 0.05 suggests a violation of the PH assumption):"
      )
      print(zph_result)

      # Save plots for Schoenfeld residuals
      message("Saving Schoenfeld residual plots to: ", output_folder)
      plot_zph <- ggcoxzph(zph_result)
      ggsave(
        str_glue("{output_folder}/ggcoxzph.pdf"),
        arrangeGrob(grobs = plot_zph)
      )

      message(
        "\nInterpretation: Look for non-random patterns in the plots (e.g., slopes, curves)."
      )
      message("A flat line around zero supports the PH assumption.")
      message("The statistical test p-value should ideally be > 0.05.")
    },
    error = function(e) {
      warning(
        "Could not perform Schoenfeld residual analysis. Error: ",
        e$message,
        call. = FALSE
      )
    }
  )

  # --- 4. Linearity Assumption (for continuous predictors) ---
  message("\n--- Checking Linearity Assumption (for continuous predictors) ---")
  message("Plotting martingale residuals against continuous predictors...")

  # --- MODIFIED SECTION ---
  tryCatch(
    {
      # Get martingale residuals
      mart_res <- residuals(fit, type = "martingale")

      # Get predictor variable names from the model terms attributes
      # This robustly gets *base variable names* (e.g., "age", "bili")
      # even if they are inside functions like ridge() or pspline()
      term_data_classes <- attr(fit$terms, "dataClasses")
      response_index <- attr(fit$terms, "response")

      # Exclude the response variable(s)
      if (!is.null(response_index) && response_index > 0) {
        predictor_vars <- names(term_data_classes)[-response_index]
      } else {
        predictor_vars <- names(term_data_classes)
      }

      if (length(predictor_vars) == 0) {
        message(
          "No predictor variables found in model terms for linearity check."
        )
        return()
      }

      # Loop through *base* predictor variables
      for (var_name in predictor_vars) {
        if (!var_name %in% names(data)) {
          message(
            "Skipping martingale plot for term '",
            var_name,
            "': variable not found in provided data."
          )
          next
        }

        var_data <- data[[var_name]]

        # Plot only if the variable is numeric (a proxy for continuous)
        if (is.numeric(var_data)) {
          # Check if it's not just a binary/categorical variable
          # We use > 5 unique values as a heuristic for "continuous"
          if (length(unique(var_data[!is.na(var_data)])) > 5) {
            message(
              "Generating martingale plot for continuous predictor: ",
              var_name
            )

            filepath <- file.path(
              output_folder,
              paste0("martingale_plot_", var_name, ".png")
            )
            png(filename = filepath, width = 800, height = 600)

            plot(
              var_data,
              mart_res,
              main = paste("Martingale Residuals vs.", var_name),
              xlab = var_name,
              ylab = "Martingale Residuals",
              pch = 19,
              col = rgb(0, 0, 0, 0.3)
            )

            # Ensure no NAs are passed to lowess
            valid_indices <- !is.na(var_data) & !is.na(mart_res)
            if (sum(valid_indices) > 2) {
              # lowess needs at least 2 points
              lines(
                lowess(var_data[valid_indices], mart_res[valid_indices]),
                col = "red",
                lwd = 2
              )
            }

            dev.off()
          } else {
            message(
              "Skipping martingale plot for '",
              var_name,
              "': appears to be categorical or binary (<= 5 unique values)."
            )
          }
        } else {
          message(
            "Skipping martingale plot for '",
            var_name,
            "': not a numeric variable."
          )
        }
      }
      # --- END MODIFIED SECTION ---

      message(
        "\nInterpretation: The red 'lowess' line should be roughly linear (a straight line)."
      )
      message(
        "Strong curves or U-shapes suggest the variable's relationship with the log-hazard is not linear,"
      )
      message(
        "and you might need to transform it (e.g., log, sqrt) or use splines."
      )
    },
    error = function(e) {
      warning(
        "Could not perform martingale residual analysis. Error: ",
        e$message,
        call. = FALSE
      )
    }
  )

  message("\n--- Assumption checks complete. ---")
  message("Plots have been saved to: ", normalizePath(output_folder))
}

# Run Elastic Net Cox regression with bootstrapping
run_enet_cox_boot <- function(
  data,
  covariates,
  time_var,
  outcome_var,
  comp_name,
  alpha = 0.5,
  n_boot = 1000,
  seed_fit = 123,
  seed_boot = 456
) {
  model_vars <- c(time_var, outcome_var, covariates)
  df_complete <- data %>%
    dplyr::select(dplyr::all_of(model_vars)) %>%
    dplyr::filter(complete.cases(.))

  x_formula <- as.formula(paste("~", paste(covariates, collapse = " + ")))
  x_mat <- model.matrix(x_formula, df_complete)[, -1]

  y_surv <- survival::Surv(df_complete[[time_var]], df_complete[[outcome_var]])

  message(sprintf("\n--- DEBUG INFO (Original Fit) for %s ---", comp_name))
  message("Number of complete cases: ", nrow(df_complete))
  message("Number of events: ", sum(df_complete[[outcome_var]]))
  message("Number of predictors in x_mat: ", ncol(x_mat))
  message("------------------------------------------")

  set.seed(seed_fit)
  withCallingHandlers(
    {
      cv_fit <- glmnet::cv.glmnet(
        x = x_mat,
        y = y_surv,
        family = "cox",
        alpha = alpha
      )
    },
    warning = function(w) {
      message(sprintf("WARNING in cv.glmnet for %s: %s", comp_name, w$message))
    }
  )

  # Save original coefficients from the non-bootstrapped dataset
  orig_coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
  orig_coefs_df <- data.frame(
    term = rownames(orig_coefs),
    estimate = orig_coefs[, 1],
    stringsAsFactors = FALSE
  )

  # Bootstrapping
  n_samples <- nrow(df_complete)
  boot_coefs <- matrix(NA, nrow = n_boot, ncol = length(orig_coefs[, 1]))
  colnames(boot_coefs) <- rownames(orig_coefs)

  set.seed(seed_boot)
  for (b in 1:n_boot) {
    boot_idx <- sample(1:n_samples, replace = TRUE)
    x_boot <- x_mat[boot_idx, , drop = FALSE]
    y_boot <- y_surv[boot_idx, ]

    tryCatch(
      {
        withCallingHandlers(
          {
            boot_fit <- glmnet::glmnet(
              x = x_boot,
              y = y_boot,
              family = "cox",
              alpha = alpha,
              maxit = 1000000
            )
            boot_coefs[b, ] <- as.matrix(coef(
              boot_fit,
              s = cv_fit$lambda.min
            ))[, 1]
          },
          warning = function(w) {
            message(sprintf(
              "Bootstrap %d WARNING: %s | Events: %d",
              b,
              w$message,
              sum(y_boot[, "status"])
            ))
          }
        )
      },
      error = function(e) {
        message(sprintf(
          "Bootstrap %d ERROR: %s | Events: %d",
          b,
          e$message,
          sum(y_boot[, "status"])
        ))
      }
    )
  }

  # Report intervals
  ci_lower <- apply(boot_coefs, 2, quantile, probs = 0.025, na.rm = TRUE)
  ci_upper <- apply(boot_coefs, 2, quantile, probs = 0.975, na.rm = TRUE)

  comp_enet_res <- orig_coefs_df %>%
    dplyr::mutate(
      dataset = comp_name,
      composite = comp_name,
      hr = exp(estimate),
      hr_low = exp(ci_lower[term]),
      hr_high = exp(ci_upper[term]),
      p.value = NA_real_,
      n_events = sum(df_complete[[outcome_var]])
    )

  return(comp_enet_res)
}


# ============================================================================
# MICE Multiple Imputation Functions
# ============================================================================

#' Impute missing data using MICE (Multiple Imputation by Chained Equations)
#'
#' @param df Data frame to impute
#' @param m Number of imputed datasets (default 20)
#' @param seed Random seed for reproducibility
#' @param exclude_from_imputation Character vector of columns to include as
#'   predictors in the imputation model but NOT impute (method = "").
#'   Typically outcome vars, auxiliary vars assumed complete.
#' @param remove_from_model Character vector of columns to exclude entirely
#'   from the MICE predictor matrix (row+col zeroed). Column is retained in
#'   completed datasets but not used during imputation. Typically ID columns.
#' @param maxit Number of MICE iterations (default 20)
#' @return A mids object from mice::mice()
mice_impute_data <- function(
  df,
  m = 20,
  seed = 123,
  exclude_from_imputation = NULL,
  remove_from_model = NULL,
  maxit = 20
) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required. Install with install.packages('mice').")
  }

  # Dry run to get default method vector and predictor matrix
  ini <- mice::mice(df, maxit = 0, printFlag = FALSE)
  method <- ini$method
  pred_matrix <- ini$predictorMatrix

  # Set method = "" for variables that should be predictors only (not imputed)
  if (!is.null(exclude_from_imputation)) {
    exclude_cols <- intersect(exclude_from_imputation, names(df))
    method[exclude_cols] <- ""
  }

  # Zero out predictor matrix rows + cols for variables to remove entirely
  if (!is.null(remove_from_model)) {
    remove_cols <- intersect(remove_from_model, names(df))
    pred_matrix[remove_cols, ] <- 0
    pred_matrix[, remove_cols] <- 0
    method[remove_cols] <- ""
  }

  # Run MICE
  mids_obj <- mice::mice(
    df,
    m = m,
    seed = seed,
    method = method,
    predictorMatrix = pred_matrix,
    maxit = maxit,
    printFlag = FALSE
  )

  return(mids_obj)
}


#' Generate MICE convergence diagnostics (trace plots + missing data summary)
#'
#' @param mids_obj A mids object from mice::mice()
#' @param output_dir Directory to save trace plot PDF
#' @param filename_prefix Prefix for the output file
#' @return invisible(mids_obj) for piping
mice_convergence_check <- function(
  mids_obj,
  output_dir,
  filename_prefix = "mice_convergence"
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Identify variables that were actually imputed (have missing data AND a method)
  has_missing <- names(which(mids_obj$nmis > 0))
  has_method <- names(which(mids_obj$method != ""))
  imputed_vars <- intersect(has_missing, has_method)

  if (length(imputed_vars) == 0) {
    message("No variables were imputed - skipping convergence check.")
    return(invisible(mids_obj))
  }

  # Print summary of missing data
  miss_summary <- data.frame(
    Variable = names(mids_obj$nmis),
    N_Missing = mids_obj$nmis,
    Pct_Missing = round(100 * mids_obj$nmis / nrow(mids_obj$data), 1),
    stringsAsFactors = FALSE
  )
  miss_summary <- miss_summary[miss_summary$N_Missing > 0, ]
  message("Missing data summary:")
  print(miss_summary, row.names = FALSE)

  # Save trace plots as PDF (handles multi-page output naturally)
  fn <- file.path(output_dir, paste0(filename_prefix, ".pdf"))
  grDevices::pdf(fn, width = 12, height = 8)
  print(plot(mids_obj))
  grDevices::dev.off()

  message(paste("Convergence trace plots saved to:", fn))
  return(invisible(mids_obj))
}


#' Pool Cox regression results across multiply imputed datasets using Rubin's Rules
#'
#' @param model_list List of m coxph model objects (from cox_function with return_model=TRUE)
#' @param predictor_name Name of the main predictor
#' @param predictor_type "continuous" or "categorical"
#' @param outcome_var Name of the event indicator variable
#' @param time_var Name of the time-to-event variable
#' @param suffix Suffix string for composite naming
#' @return Tibble matching cox_function() output format with pooled estimates
pool_cox_rubin <- function(
  model_list,
  predictor_name,
  predictor_type,
  outcome_var,
  time_var,
  suffix,
  entry_var = NULL # entry-time variable name (left-truncation); metadata only
) {
  m <- length(model_list)
  stopifnot(m >= 2)

  # Extract tidy results from each model on log-HR scale
  tidy_list <- lapply(model_list, function(fit) {
    broom::tidy(fit, exponentiate = FALSE) |>
      dplyr::select(term, estimate, std.error)
  })

  # Collect all unique terms across imputations
  all_terms <- unique(unlist(lapply(tidy_list, function(x) x$term)))

  # Pool each term via Rubin's Rules
  pooled_rows <- lapply(all_terms, function(t) {
    betas <- numeric(m)
    ses <- numeric(m)
    valid <- logical(m)

    for (i in seq_len(m)) {
      row <- tidy_list[[i]] |> dplyr::filter(term == t)
      if (nrow(row) == 1) {
        betas[i] <- row$estimate
        ses[i] <- row$std.error
        valid[i] <- TRUE
      }
    }

    if (sum(valid) < 2) {
      return(tibble::tibble(
        term = t,
        hr = NA_real_,
        hr_low = NA_real_,
        hr_high = NA_real_,
        p.value = NA_real_
      ))
    }

    betas_v <- betas[valid]
    ses_v <- ses[valid]
    m_eff <- length(betas_v)

    pooled_beta <- mean(betas_v)
    W <- mean(ses_v^2) # within-imputation variance
    B <- stats::var(betas_v) # between-imputation variance
    T_var <- W + (1 + 1 / m_eff) * B # total variance
    pooled_se <- sqrt(T_var)

    # Barnard-Rubin degrees of freedom
    if (B < .Machine$double.eps) {
      df_val <- Inf
    } else {
      lambda <- (1 + 1 / m_eff) * B / T_var
      df_val <- (m_eff - 1) / lambda^2
    }

    t_stat <- pooled_beta / pooled_se
    p_val <- 2 * stats::pt(abs(t_stat), df = df_val, lower.tail = FALSE)
    t_crit <- stats::qt(0.975, df = df_val)

    tibble::tibble(
      term = t,
      hr = exp(pooled_beta),
      hr_low = exp(pooled_beta - t_crit * pooled_se),
      hr_high = exp(pooled_beta + t_crit * pooled_se),
      p.value = p_val
    )
  })

  pooled_tbl <- dplyr::bind_rows(pooled_rows)

  # Extract metadata from first model
  first_fit <- model_list[[1]]
  ph_failed_vars <- attr(first_fit, "ph_failed_vars")
  if (is.null(ph_failed_vars)) {
    ph_failed_vars <- character(0)
  }

  td_str <- if (length(ph_failed_vars) == 0) {
    NA_character_
  } else {
    paste(ph_failed_vars, collapse = ",")
  }

  out <- pooled_tbl |>
    dplyr::mutate(
      predictor = predictor_name,
      predictor_type = predictor_type,
      time_dependent_predictor = td_str,
      n = stats::nobs(first_fit),
      n_events = first_fit$nevent,
      outcome_entry_var = if (is.null(entry_var)) NA_character_ else entry_var,
      outcome_time_var = time_var,
      outcome_event_var = outcome_var,
      composite = stringr::str_remove(suffix, "^_"),
      .before = 1
    )

  attr(out, "ph_failed_vars") <- ph_failed_vars
  out
}


#' Run a fixed multivariable Cox model across multiply imputed datasets and pool via Rubin's Rules
#'
#' Operates on a mids object. A single, fixed (plain) Cox model with an identical
#' formula is fit on each of the m completed datasets and the results are pooled
#' via Rubin's Rules (pool_cox_rubin). No PH check or time-dependent (tt()) terms
#' are applied: an adaptive specification could differ across imputations, which
#' would break pooling. Assess PH separately as a per-imputation diagnostic
#' (see ph_diagnostic_mice).
#'
#' @param mids_obj A mids object from mice_impute_data()
#' @param predictor Main predictor variable name
#' @param predictor_type "continuous" or "categorical"
#' @param covariates Character vector of covariate names
#' @param outcome Event indicator variable name
#' @param time Time-to-event variable name
#' @param suffix Suffix for output naming
#' @return Pooled tibble matching run_adaptive_cox() output format

run_adaptive_cox_mice <- function(
  mids_obj,
  predictor,
  predictor_type,
  covariates,
  outcome,
  time,
  suffix,
  entry_time = NULL # entry-time variable for left-truncation/delayed entry; passed through to cox_function
) {
  m <- mids_obj$m

  # Extract completed datasets
  completed_datasets <- lapply(seq_len(m), function(i) {
    mice::complete(mids_obj, i)
  })

  # Fit a fixed plain Cox model (identical formula) on each imputed dataset.
  model_list <- lapply(completed_datasets, function(d) {
    cox_function(
      input_data = d,
      predictor_name = predictor,
      predictor_type = predictor_type,
      covars = covariates,
      indicator_varname = outcome,
      time_varname = time,
      entry_time_varname = entry_time,
      output_dir = NULL,
      suffix = suffix,
      return_model = TRUE,
      run_PHassumption_check = FALSE
    )
  })

  # Pool via Rubin's Rules
  pool_cox_rubin(
    model_list = model_list,
    predictor_name = predictor,
    predictor_type = predictor_type,
    outcome_var = outcome,
    time_var = time,
    entry_var = entry_time,
    suffix = suffix
  )
}


#' Proportional-hazards diagnostic across multiply imputed datasets
#'
#' Fits a fixed (plain) Cox model with an identical formula on each of the m
#' completed datasets, runs cox.zph on each, and reports how many of the m
#' imputations flag a PH violation (p < 0.05) per term. cox.zph cannot be
#' Rubin-pooled cleanly, so this is reported as a diagnostic only; the model
#' specification is kept fixed regardless of the result.
#'
#' @param mids_obj A mids object
#' @param predictor Main predictor variable name
#' @param covariates Character vector of covariate names
#' @param outcome Event indicator variable name
#' @param time Time-to-event variable name
#' @param alpha Per-imputation significance threshold for flagging (default 0.05)
#' @return Tibble: term, m, n_flagged, pct_flagged (descending by n_flagged)
ph_diagnostic_mice <- function(
  mids_obj,
  predictor,
  covariates,
  outcome,
  time,
  entry_time = NULL, # entry-time variable for left-truncation; mirror the fitted model's Surv
  alpha = 0.05
) {
  m <- mids_obj$m
  completed <- lapply(seq_len(m), function(i) mice::complete(mids_obj, i))

  rhs <- paste(unique(c(predictor, covariates)), collapse = " + ")
  surv_lhs <- if (is.null(entry_time)) {
    paste0("Surv(", time, ", ", outcome, ")")
  } else {
    paste0("Surv(", entry_time, ", ", time, ", ", outcome, ")")
  }
  fml <- stats::as.formula(paste0(surv_lhs, " ~ ", rhs))

  zph_list <- lapply(completed, function(d) {
    out <- tryCatch(
      {
        fit <- survival::coxph(fml, data = d)
        z <- survival::cox.zph(fit)
        tab <- as.data.frame(z$table)
        tab$term <- rownames(tab)
        tab[tab$term != "GLOBAL", c("term", "p"), drop = FALSE]
      },
      error = function(e) NULL
    )
    out
  })

  zph_all <- dplyr::bind_rows(zph_list)
  if (nrow(zph_all) == 0) {
    return(tibble::tibble(
      term = character(0),
      m = integer(0),
      n_flagged = integer(0),
      pct_flagged = numeric(0)
    ))
  }

  zph_all |>
    dplyr::group_by(term) |>
    dplyr::summarise(
      m = dplyr::n(),
      n_flagged = sum(p < alpha, na.rm = TRUE),
      pct_flagged = round(100 * sum(p < alpha, na.rm = TRUE) / dplyr::n(), 1),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_flagged))
}


#' Elastic Net Cox regression with bootstrapping across multiply imputed datasets
#'
#' Uses median lambda strategy: find lambda.min via CV in each imputed dataset,
#' take the median, then fit + bootstrap at that lambda across all imputations.
#'
#' @param mids_obj A mids object from mice_impute_data()
#' @param covariates Character vector of covariate/predictor names
#' @param time_var Time-to-event variable name
#' @param outcome_var Event indicator variable name
#' @param comp_name Name of the composite/dataset
#' @param alpha Elastic net mixing parameter (default 0.5)
#' @param n_boot Total number of bootstrap replicates (split across imputations)
#' @param seed_fit Seed for CV fitting
#' @param seed_boot Seed for bootstrapping
#' @return Tibble matching run_enet_cox_boot() output format
run_enet_cox_boot_mice <- function(
  mids_obj,
  covariates,
  time_var,
  outcome_var,
  comp_name,
  alpha = 0.5,
  n_boot = 1000,
  seed_fit = 123,
  seed_boot = 456
) {
  m <- mids_obj$m
  completed_datasets <- lapply(1:m, function(i) mice::complete(mids_obj, i))

  model_vars <- c(time_var, outcome_var, covariates)

  # Prepare model matrices for each imputed dataset
  prepared <- lapply(completed_datasets, function(df) {
    df_complete <- df |>
      dplyr::select(dplyr::all_of(model_vars)) |>
      dplyr::filter(stats::complete.cases(dplyr::across(dplyr::everything())))

    x_formula <- stats::as.formula(paste(
      "~",
      paste(covariates, collapse = " + ")
    ))
    x_mat <- stats::model.matrix(x_formula, df_complete)[, -1]
    y_surv <- survival::Surv(
      df_complete[[time_var]],
      df_complete[[outcome_var]]
    )

    list(x_mat = x_mat, y_surv = y_surv, df_complete = df_complete)
  })

  # Step 1: CV for lambda in each imputed dataset
  set.seed(seed_fit)
  lambda_mins <- numeric(m)

  for (i in seq_len(m)) {
    withCallingHandlers(
      {
        cv_fit <- glmnet::cv.glmnet(
          x = prepared[[i]]$x_mat,
          y = prepared[[i]]$y_surv,
          family = "cox",
          alpha = alpha
        )
        lambda_mins[i] <- cv_fit$lambda.min
      },
      warning = function(w) {
        message(sprintf(
          "WARNING in cv.glmnet for %s (imputation %d): %s",
          comp_name,
          i,
          w$message
        ))
      }
    )
  }

  median_lambda <- stats::median(lambda_mins)
  message(sprintf(
    "Median lambda across %d imputations for %s: %.6f",
    m,
    comp_name,
    median_lambda
  ))

  # Step 2: Fit at median lambda on each imputed dataset
  coef_list <- lapply(seq_len(m), function(i) {
    fit <- glmnet::glmnet(
      x = prepared[[i]]$x_mat,
      y = prepared[[i]]$y_surv,
      family = "cox",
      alpha = alpha
    )
    as.matrix(stats::coef(fit, s = median_lambda))[, 1]
  })

  coef_mat <- do.call(rbind, coef_list)
  avg_coefs <- colMeans(coef_mat)

  # Step 3: Bootstrap across imputed datasets
  n_boot_per_imp <- ceiling(n_boot / m)
  total_boot <- n_boot_per_imp * m

  boot_coefs <- matrix(NA, nrow = total_boot, ncol = length(avg_coefs))
  colnames(boot_coefs) <- names(avg_coefs)

  set.seed(seed_boot)
  boot_idx <- 0

  for (i in seq_len(m)) {
    x_mat <- prepared[[i]]$x_mat
    y_surv <- prepared[[i]]$y_surv
    n_samples <- nrow(x_mat)

    for (b in seq_len(n_boot_per_imp)) {
      boot_idx <- boot_idx + 1
      idx <- sample(seq_len(n_samples), replace = TRUE)
      x_boot <- x_mat[idx, , drop = FALSE]
      y_boot <- y_surv[idx, ]

      tryCatch(
        {
          suppressWarnings({
            boot_fit <- glmnet::glmnet(
              x = x_boot,
              y = y_boot,
              family = "cox",
              alpha = alpha,
              maxit = 1000000
            )
            boot_coefs[boot_idx, ] <- as.matrix(stats::coef(
              boot_fit,
              s = median_lambda
            ))[, 1]
          })
        },
        error = function(e) {
          message(sprintf(
            "Bootstrap %d (imp %d) ERROR for %s: %s",
            b,
            i,
            comp_name,
            e$message
          ))
        }
      )
    }
  }

  # Step 4: Compute CIs from pooled bootstrap distribution
  ci_lower <- apply(boot_coefs, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  ci_upper <- apply(boot_coefs, 2, stats::quantile, probs = 0.975, na.rm = TRUE)

  # Step 5: Format output to match run_enet_cox_boot()
  n_events_first <- sum(prepared[[1]]$df_complete[[outcome_var]])

  result <- data.frame(
    term = names(avg_coefs),
    estimate = avg_coefs,
    stringsAsFactors = FALSE
  ) |>
    dplyr::mutate(
      dataset = comp_name,
      composite = comp_name,
      hr = exp(estimate),
      hr_low = exp(ci_lower[term]),
      hr_high = exp(ci_upper[term]),
      p.value = NA_real_,
      n_events = n_events_first
    )

  return(result)
}

# Events Per Variable
# EPV = Number of Events / Number of Predictors
# Predictors = Main Predictor + Covariates (degrees of freedom approx)

calculate_epv <- function(model_tibble, n_predictors) {
  model_tibble %>%
    select(dataset, term, n_events, n) %>%
    distinct(dataset, n_events, n) %>%
    mutate(
      n_predictors = n_predictors,
      epv = n_events / n_predictors
    )
}
