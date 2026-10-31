# ===== PLOTTING HELPERS =====
#
# Functions for generating QQ and volcano plots from FAST_outcome_WAS results.

library(qqman)
library(ggplot2)

.coerce_to_strata_results <- function(results, analysis = c("analysis_change", "analysis_level")) {
  analysis <- match.arg(analysis)

  if (is.list(results) && all(c("all", "male", "female") %in% names(results))) {
    return(results)
  }

  if (is.list(results) && analysis %in% names(results)) {
    return(results[[analysis]])
  }

  stop("Invalid results object. Pass FAST_outcome_WAS output or a stratum-level results list.")
}


.outcome_filtered_membership <- function(outcome_effects) {
  if ("IN_FILTERED_SET" %in% names(outcome_effects)) {
    return(outcome_effects$IN_FILTERED_SET %in% TRUE)
  }
  # Older results encoded membership only through available filtered inference.
  if ("BH_P_VALUE_FILTERED" %in% names(outcome_effects)) {
    p <- outcome_effects$BH_P_VALUE_FILTERED
    return(is.finite(p) & p >= 0 & p <= 1)
  }
  rep(FALSE, nrow(outcome_effects))
}


.outcome_plot_rows <- function(results, stratum, fu, probe_set) {
  oe <- results[[stratum]]$outcome_effects
  if (is.null(oe)) stop(paste0("No results for stratum '", stratum, "'"))
  if ("ARM" %in% names(oe)) {
    stop("Outcome plots require one common analyte slope per ANALYTE_NAME x FU; ARM tables are not supported.")
  }
  oe <- oe[oe$FU == fu, , drop = FALSE]
  if (nrow(oe) == 0L) stop(paste0("No results for FU=", fu))
  if (probe_set == "filtered") {
    oe <- oe[.outcome_filtered_membership(oe), , drop = FALSE]
  }
  oe
}


.outcome_plot_estimable <- function(oe) {
  keep <- is.finite(oe$P_VALUE) & oe$P_VALUE >= 0 & oe$P_VALUE <= 1 &
    is.finite(oe$EFFECT_SIZE)
  if ("STATUS" %in% names(oe)) keep <- keep & oe$STATUS %in% "ok"
  keep
}

.outcome_fu_plot_label <- function(outcome_effects, fu) {
  if (!"FU_LABEL" %in% names(outcome_effects)) return(paste0("FU=", fu))
  labels <- unique(as.character(outcome_effects$FU_LABEL[outcome_effects$FU == fu]))
  if (length(labels) != 1L || anyNA(labels) || !nzchar(trimws(labels))) {
    stop("Each FU must identify exactly one non-blank FU_LABEL")
  }
  labels
}


.outcome_empty_plot_label <- function(oe, probe_set) {
  if (nrow(oe) == 0L && probe_set == "filtered") {
    "No probes in filtered set"
  } else {
    "No estimable results"
  }
}


plot_qq <- function(results, stratum = "all", fu = 1, probe_set = "full",
                    analysis = c("analysis_change", "analysis_level"), title = NULL) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))
  probe_set <- match.arg(probe_set, c("full", "filtered"))
  oe <- .outcome_plot_rows(results, stratum, fu, probe_set)
  pvals <- oe$P_VALUE[.outcome_plot_estimable(oe)]

  if (is.null(title)) {
    title <- paste0("QQ Plot - ", stratum, ", ",
                    .outcome_fu_plot_label(results[[stratum]]$outcome_effects, fu), ", ", probe_set)
  }

  if (!length(pvals)) {
    plot.new()
    title(main = title)
    text(0.5, 0.5, .outcome_empty_plot_label(oe, probe_set))
    return(invisible(NULL))
  }
  qq(pmax(pvals, .Machine$double.xmin), main = title)
}


plot_volcano <- function(results, stratum = "all", fu = 1, probe_set = "full",
                         p_threshold = 0.05,
                         analysis = c("analysis_change", "analysis_level"),
                         title = NULL) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))
  probe_set <- match.arg(probe_set, c("full", "filtered"))
  oe <- .outcome_plot_rows(results, stratum, fu, probe_set)
  if (is.null(title)) {
    title <- paste0("Volcano Plot - ", stratum, ", ",
                    .outcome_fu_plot_label(results[[stratum]]$outcome_effects, fu), ", ", probe_set)
  }
  keep <- .outcome_plot_estimable(oe)
  if (!any(keep)) {
    return(ggplot() +
      annotate("text", x = 0, y = 0, label = .outcome_empty_plot_label(oe, probe_set)) +
      labs(title = title) + theme_void())
  }
  oe <- oe[keep, , drop = FALSE]

  bh_col <- if (probe_set == "filtered") "BH_P_VALUE_FILTERED" else "BH_P_VALUE"
  sig <- !is.na(oe[[bh_col]]) & oe[[bh_col]] < p_threshold

  oe$neg_log_p <- -log10(pmax(oe$P_VALUE, .Machine$double.xmin))
  oe$direction <- ifelse(!sig, "Not significant",
                         ifelse(oe$EFFECT_SIZE > 0, "Up", "Down"))

  all_levels <- c("Down", "Not significant", "Up")
  oe$direction <- factor(oe$direction, levels = all_levels)

  for (lvl in setdiff(all_levels, as.character(unique(oe$direction)))) {
    phantom <- oe[NA, ]
    phantom$direction <- factor(lvl, levels = all_levels)
    oe <- rbind(oe, phantom)
  }

  ggplot(oe, aes(x = EFFECT_SIZE, y = neg_log_p, color = direction)) +
    geom_point(size = 1, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = c("Down" = "blue", "Not significant" = "grey50", "Up" = "red")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    labs(x = "Effect Size", y = expression(-log[10](p)), title = title, color = NULL) +
    theme_minimal() +
    theme(legend.position = "right")
}


generate_all_plots <- function(results, figures_dir = NULL,
                               analysis = c("analysis_change", "analysis_level")) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  if (is.null(figures_dir)) {
    figures_dir <- "Figures"
    if (!dir.exists(figures_dir)) {
      response <- readline(paste0("'", figures_dir, "' does not exist. Create it? (y/n): "))
      if (!tolower(response) %in% c("y", "yes")) {
        stop("No output directory. Provide figures_dir or create 'Figures/'.")
      }
      dir.create(figures_dir)
    }
  } else if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
  }

  strata <- c("all", "male", "female")[vapply(results[c("all", "male", "female")],
    function(x) !is.null(x$outcome_effects) && nrow(x$outcome_effects) > 0L, logical(1))]
  if (!length(strata)) {
    message("No outcome-effect rows to plot.")
    return(invisible(NULL))
  }
  fu_by_stratum <- setNames(lapply(strata, function(s) {
    sort(unique(results[[s]]$outcome_effects$FU))
  }), strata)
  fu_levels <- sort(unique(unlist(fu_by_stratum, use.names = FALSE)))
  has_filtered <- any(vapply(strata, function(s) {
    oe <- results[[s]]$outcome_effects
    "IN_FILTERED_SET" %in% names(oe) || any(.outcome_filtered_membership(oe))
  }, logical(1)))
  probe_sets <- if (has_filtered) c("full", "filtered") else c("full")

  cat("Generating plots:\n")
  cat("  Strata:     ", paste(strata, collapse = ", "), "\n")
  cat("  FU levels:  ", paste(fu_levels, collapse = ", "), "\n")
  cat("  Probe sets: ", paste(probe_sets, collapse = ", "), "\n")

  save_plot_pdf <- function(path, width, height, draw) {
    pdf(path, width = width, height = height)
    on.exit(dev.off(), add = TRUE)
    draw()
  }

  for (ps in probe_sets) {
    qq_file <- file.path(figures_dir, paste0("qq_", ps, ".pdf"))
    save_plot_pdf(qq_file, 7, 6, function() {
      for (s in strata) {
        for (fu in fu_by_stratum[[s]]) {
          plot_qq(results, stratum = s, fu = fu, probe_set = ps)
        }
      }
    })
    cat("  Saved:", qq_file, "\n")

    volcano_file <- file.path(figures_dir, paste0("volcano_", ps, ".pdf"))
    save_plot_pdf(volcano_file, 8, 6, function() {
      for (s in strata) {
        for (fu in fu_by_stratum[[s]]) {
          print(plot_volcano(results, stratum = s, fu = fu, probe_set = ps))
        }
      }
    })
    cat("  Saved:", volcano_file, "\n")
  }

  cat("Done.\n")
}
