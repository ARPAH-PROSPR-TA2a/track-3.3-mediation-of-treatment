.adjust_within_fu <- function(p_values, fu) {
  adjusted <- rep(NA_real_, length(p_values))
  for (level in sort(unique(fu))) {
    index <- which(fu == level & !is.na(p_values))
    adjusted[index] <- p.adjust(p_values[index], method = "BH")
  }
  adjusted
}


.compute_mediation <- function(treatment_entry, outcome_entry, analysis, stratum,
                               outcome_type) {
  treatment <- treatment_entry$treatment_effects
  outcome <- outcome_entry$outcome_effects

  treatment_filtered <- "BH_P_VALUE_FILTERED" %in% names(treatment)
  outcome_filtered <- "BH_P_VALUE_FILTERED" %in% names(outcome)
  if (xor(treatment_filtered, outcome_filtered)) {
    stop(
      "Filtered-probe columns disagree for ", analysis, "/", stratum,
      ": both inputs must either contain BH_P_VALUE_FILTERED or omit it."
    )
  }

  alpha <- data.frame(
    ANALYTE_NAME = as.character(treatment$ANALYTE_NAME),
    FU = as.integer(as.character(treatment$FU)),
    ALPHA_EFFECT = treatment$EFFECT_SIZE,
    ALPHA_SE = treatment$SE,
    stringsAsFactors = FALSE
  )
  beta <- data.frame(
    ANALYTE_NAME = as.character(outcome$ANALYTE_NAME),
    FU = as.integer(as.character(outcome$FU)),
    BETA_EFFECT = outcome$EFFECT_SIZE,
    BETA_SE = outcome$SE,
    stringsAsFactors = FALSE
  )

  if (treatment_filtered) {
    alpha$ALPHA_FILTERED <- !is.na(treatment$BH_P_VALUE_FILTERED)
    beta$BETA_FILTERED <- !is.na(outcome$BH_P_VALUE_FILTERED)
  }

  alpha_key <- paste(alpha$ANALYTE_NAME, alpha$FU, sep = "\r")
  beta_key <- paste(beta$ANALYTE_NAME, beta$FU, sep = "\r")
  shared_key <- intersect(alpha_key, beta_key)
  diagnostics <- list(
    treatment_rows = nrow(alpha),
    outcome_rows = nrow(beta),
    matched_rows = length(shared_key),
    treatment_only_rows = sum(!alpha_key %in% beta_key),
    outcome_only_rows = sum(!beta_key %in% alpha_key)
  )

  if (length(shared_key) == 0L) {
    stop("No matching ANALYTE_NAME x FU keys for ", analysis, "/", stratum, ".")
  }
  if (diagnostics$treatment_only_rows > 0L || diagnostics$outcome_only_rows > 0L) {
    warning(
      "Partial mediation join for ", analysis, "/", stratum, ": ",
      diagnostics$matched_rows, " matched, ",
      diagnostics$treatment_only_rows, " treatment-only, and ",
      diagnostics$outcome_only_rows, " outcome-only rows."
    )
  }

  mediation <- merge(alpha, beta, by = c("ANALYTE_NAME", "FU"), sort = FALSE)

  if (treatment_filtered && any(mediation$ALPHA_FILTERED != mediation$BETA_FILTERED)) {
    stop("Filtered-probe membership disagrees for ", analysis, "/", stratum, ".")
  }

  mediation$INDIRECT_EFFECT <- mediation$ALPHA_EFFECT * mediation$BETA_EFFECT
  mediation$INDIRECT_SE <- sqrt(
    mediation$BETA_EFFECT^2 * mediation$ALPHA_SE^2 +
      mediation$ALPHA_EFFECT^2 * mediation$BETA_SE^2
  )

  singular <- mediation$INDIRECT_SE == 0 & mediation$INDIRECT_EFFECT == 0
  mediation$Z_VALUE <- mediation$INDIRECT_EFFECT / mediation$INDIRECT_SE
  mediation$Z_VALUE[singular] <- 0
  mediation$P_VALUE <- 2 * pnorm(-abs(mediation$Z_VALUE))

  critical_value <- qnorm(0.975)
  mediation$CI_LOWER <- mediation$INDIRECT_EFFECT - critical_value * mediation$INDIRECT_SE
  mediation$CI_UPPER <- mediation$INDIRECT_EFFECT + critical_value * mediation$INDIRECT_SE
  mediation$BH_P_VALUE <- .adjust_within_fu(mediation$P_VALUE, mediation$FU)

  if (treatment_filtered) {
    filtered <- mediation$ALPHA_FILTERED
    mediation$BH_P_VALUE_FILTERED <- NA_real_
    for (level in sort(unique(mediation$FU))) {
      index <- which(mediation$FU == level & filtered)
      mediation$BH_P_VALUE_FILTERED[index] <- p.adjust(
        mediation$P_VALUE[index], method = "BH"
      )
    }
    mediation$ALPHA_FILTERED <- NULL
    mediation$BETA_FILTERED <- NULL
  }

  if (outcome_type == "tte") {
    mediation$INDIRECT_HAZARD_RATIO <- exp(mediation$INDIRECT_EFFECT)
  }

  mediation <- mediation[order(mediation$ANALYTE_NAME, mediation$FU), , drop = FALSE]
  row.names(mediation) <- NULL

  list(
    mediation_effects = mediation,
    join_diagnostics = diagnostics
  )
}
