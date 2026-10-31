source(file.path("mediation", "R", "validation_helpers.R"))
source(file.path("mediation", "R", "analysis_helpers.R"))


FAST_mediation <- function(treatment_results,
                           outcome_results,
                           outcome_type = c("continuous", "tte")) {
  outcome_type <- match.arg(outcome_type)
  treatment_results <- .validate_results_object(treatment_results, role = "treatment")
  outcome_results <- .validate_results_object(
    outcome_results,
    role = "outcome",
    outcome_type = outcome_type
  )
  fu_labels <- .resolve_mediation_fu_labels(treatment_results, outcome_results)

  outcome_labels <- .mediation_metadata_labels(outcome_results, "outcome")
  treatment_labels <- .mediation_metadata_labels(treatment_results, "treatment")
  treatment_map <- attr(treatment_results, "treatment_map")
  if (!is.null(outcome_labels)) {
    requested <- if (!is.null(treatment_labels)) names(treatment_labels) else c("0", treatment_map$ARM)
    if (any(!requested %in% names(outcome_labels))) {
      stop("Outcome treatment metadata does not include every requested alpha arm and control.")
    }
    if (!is.null(treatment_labels) && !setequal(names(treatment_labels), names(outcome_labels))) {
      stop("Alpha and beta treatment metadata describe different arm sets.")
    }
  }

  output <- setNames(vector("list", length(.supported_analyses)), .supported_analyses)

  for (analysis in .supported_analyses) {
    treatment_strata <- intersect(.supported_strata, names(treatment_results[[analysis]]))
    outcome_strata <- intersect(.supported_strata, names(outcome_results[[analysis]]))
    strata <- .supported_strata[
      .supported_strata %in% union(treatment_strata, outcome_strata)
    ]
    output[[analysis]] <- setNames(vector("list", length(strata)), strata)

    for (stratum in strata) {
      treatment_entry <- treatment_results[[analysis]][[stratum]]
      outcome_entry <- outcome_results[[analysis]][[stratum]]

      if (is.null(treatment_entry) && is.null(outcome_entry)) {
        output[[analysis]][[stratum]] <- NULL
        next
      }
      if (is.null(treatment_entry) || is.null(outcome_entry)) {
        warning(
          "Skipping ", analysis, "/", stratum,
          " because the stratum is available in only one input."
        )
        output[[analysis]][[stratum]] <- NULL
        next
      }

      output[[analysis]][[stratum]] <- .compute_mediation(
        treatment_entry = treatment_entry,
        outcome_entry = outcome_entry,
        analysis = analysis,
        stratum = stratum,
        outcome_type = outcome_type,
        outcome_labels = outcome_labels
      )
      # Presentation only: unknown visits in legacy partial inputs remain NA.
      if (!is.null(fu_labels)) {
        effects <- output[[analysis]][[stratum]]$mediation_effects
        effects$FU_LABEL <- unname(fu_labels[as.character(effects$FU)])
        columns <- setdiff(names(effects), "FU_LABEL")
        output[[analysis]][[stratum]]$mediation_effects <- effects[
          , append(columns, "FU_LABEL", after = match("FU", columns)), drop = FALSE
        ]
      }
    }
  }

  attr(output, "outcome_type") <- outcome_type
  attr(output, "method") <- "product_of_coefficients"
  attr(output, "inference") <- "sobel"
  if (is.null(treatment_labels) && !is.null(treatment_map)) {
    treatment_labels <- c(`0` = treatment_map$REFERENCE[1L],
                          setNames(treatment_map$TREATMENT, treatment_map$ARM))
  }
  attr(output, "metadata") <- list(
    schema_version = "2.0", analyte_slope = "common_across_treatment_arms",
    treatment_labels = treatment_labels, treatments = treatment_map,
    estimand = if (outcome_type == "tte") "cox_product_of_coefficients_approximation" else {
      "pure_and_total_natural_indirect_effect_under_shared_beta"
    },
    effect_scale = if (outcome_type == "tte") "log_hazard" else "outcome",
    bh_family = "BH within analysis, stratum, FU and contrast",
    alpha_beta_covariance = "not_included", cross_contrast_covariance = "not_available",
    survival_interpretation = if (outcome_type == "tte") {
      "exp(product) is a Cox product approximation, not a general causal indirect hazard ratio"
    } else NULL
  )
  if (!is.null(fu_labels)) attr(output, "metadata")$fu_labels <- fu_labels
  output
}
