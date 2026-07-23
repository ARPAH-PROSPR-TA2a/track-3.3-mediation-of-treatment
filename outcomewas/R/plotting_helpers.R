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


plot_qq <- function(results, stratum = "all", fu = 1, probe_set = "full",
                    analysis = c("analysis_change", "analysis_level"), title = NULL) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  oe <- results[[stratum]]$outcome_effects
  if (is.null(oe)) stop(paste0("No results for stratum '", stratum, "'"))

  oe <- oe[oe$FU == fu, ]
  if (nrow(oe) == 0) stop(paste0("No results for FU=", fu))

  if (probe_set == "filtered") {
    oe <- oe[!is.na(oe$BH_P_VALUE_FILTERED), ]
    if (nrow(oe) == 0) stop("No filtered probes found")
  }

  pvals <- oe$P_VALUE[!is.na(oe$P_VALUE)]

  if (is.null(title)) {
    title <- paste0("QQ Plot - ", stratum, ", FU=", fu, ", ", probe_set)
  }

  qq(pvals, main = title)
}


plot_volcano <- function(results, stratum = "all", fu = 1, probe_set = "full",
                         p_threshold = 0.05,
                         analysis = c("analysis_change", "analysis_level"),
                         title = NULL) {

  results <- .coerce_to_strata_results(results, analysis = match.arg(analysis))

  oe <- results[[stratum]]$outcome_effects
  if (is.null(oe)) stop(paste0("No results for stratum '", stratum, "'"))
  oe <- oe[oe$FU == fu, ]
  if (nrow(oe) == 0) stop(paste0("No results for FU=", fu))

  if (probe_set == "filtered") {
    oe <- oe[!is.na(oe$BH_P_VALUE_FILTERED), ]
    if (nrow(oe) == 0) stop("No filtered probes found")
  }

  bh_col <- if (probe_set == "filtered") "BH_P_VALUE_FILTERED" else "BH_P_VALUE"
  sig <- !is.na(oe[[bh_col]]) & oe[[bh_col]] < p_threshold

  oe$neg_log_p <- -log10(oe$P_VALUE)
  oe$direction <- ifelse(!sig, "Not significant",
                         ifelse(oe$EFFECT_SIZE > 0, "Up", "Down"))

  all_levels <- c("Down", "Not significant", "Up")
  oe$direction <- factor(oe$direction, levels = all_levels)

  for (lvl in setdiff(all_levels, as.character(unique(oe$direction)))) {
    phantom <- oe[NA, ]
    phantom$direction <- factor(lvl, levels = all_levels)
    oe <- rbind(oe, phantom)
  }

  if (is.null(title)) {
    title <- paste0("Volcano Plot - ", stratum, ", FU=", fu, ", ", probe_set)
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

  strata <- c("all", "male", "female")[!sapply(results[c("all", "male", "female")], is.null)]
  fu_levels <- sort(unique(results[[strata[1]]]$outcome_effects$FU))

  has_filtered <- "BH_P_VALUE_FILTERED" %in% colnames(results[[strata[1]]]$outcome_effects) &&
    any(!is.na(results[[strata[1]]]$outcome_effects$BH_P_VALUE_FILTERED))
  probe_sets <- if (has_filtered) c("full", "filtered") else c("full")

  cat("Generating plots:\n")
  cat("  Strata:     ", paste(strata, collapse = ", "), "\n")
  cat("  FU levels:  ", paste(fu_levels, collapse = ", "), "\n")
  cat("  Probe sets: ", paste(probe_sets, collapse = ", "), "\n")

  for (ps in probe_sets) {
    qq_file <- file.path(figures_dir, paste0("qq_", ps, ".pdf"))
    pdf(qq_file, width = 7, height = 6)
    for (s in strata) {
      for (fu in fu_levels) {
        plot_qq(results, stratum = s, fu = fu, probe_set = ps)
      }
    }
    dev.off()
    cat("  Saved:", qq_file, "\n")

    volcano_file <- file.path(figures_dir, paste0("volcano_", ps, ".pdf"))
    pdf(volcano_file, width = 8, height = 6)
    for (s in strata) {
      for (fu in fu_levels) {
        print(plot_volcano(results, stratum = s, fu = fu, probe_set = ps))
      }
    }
    dev.off()
    cat("  Saved:", volcano_file, "\n")
  }

  cat("Done.\n")
}
