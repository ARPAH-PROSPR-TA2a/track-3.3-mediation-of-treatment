source(file.path("mediation", "R", "validation_helpers.R"))
source(file.path("mediation", "R", "analysis_helpers.R"))


FAST_mediation <- function(treatment_results,
                           outcome_results,
                           outcome_type = c("continuous", "tte")) {
  outcome_type <- match.arg(outcome_type)
  .validate_results_object(treatment_results, role = "treatment")
  .validate_results_object(
    outcome_results,
    role = "outcome",
    outcome_type = outcome_type
  )

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
        outcome_type = outcome_type
      )
    }
  }

  attr(output, "outcome_type") <- outcome_type
  attr(output, "method") <- "product_of_coefficients"
  attr(output, "inference") <- "sobel"
  output
}
