.adjust_within_fu <- function(p_values, fu, contrast = rep("1_vs_0", length(fu))) {
  adjusted <- rep(NA_real_, length(p_values))
  groups <- paste(fu, contrast, sep = "\r")
  for (group in unique(groups)) {
    index <- which(groups == group & !is.na(p_values))
    adjusted[index] <- p.adjust(p_values[index], method = "BH")
  }
  adjusted
}

.compute_mediation <- function(treatment_entry, outcome_entry, analysis, stratum,
                               outcome_type, outcome_labels = NULL) {
  treatment <- treatment_entry$treatment_effects
  outcome <- outcome_entry$outcome_effects
  treatment_filtered <- "IN_FILTERED_SET" %in% names(treatment)
  outcome_filtered <- "IN_FILTERED_SET" %in% names(outcome)
  if (xor(treatment_filtered, outcome_filtered)) {
    stop("Filtered-probe columns disagree for ", analysis, "/", stratum,
         ": both inputs must provide IN_FILTERED_SET or complete legacy BH_P_VALUE_FILTERED membership.")
  }
  alpha <- data.frame(
    ANALYTE_NAME = treatment$ANALYTE_NAME, FU = treatment$FU,
    CONTRAST = treatment$CONTRAST, TREATMENT = treatment$TREATMENT,
    REFERENCE = treatment$REFERENCE, ALPHA_EFFECT = treatment$EFFECT_SIZE,
    ALPHA_SE = treatment$SE, ALPHA_STATUS = treatment$STATUS,
    stringsAsFactors = FALSE
  )
  alpha_arm <- sub("_vs_0$", "", alpha$CONTRAST)
  alpha_key <- paste(alpha$ANALYTE_NAME, alpha$FU, sep = "\r")
  beta_key <- paste(outcome$ANALYTE_NAME, outcome$FU, sep = "\r")
  if (!is.null(outcome_labels)) {
    supplied <- attr(treatment, "supplied_labels")
    known_active <- alpha_arm %in% names(outcome_labels)
    if (("TREATMENT" %in% supplied &&
         any(alpha$TREATMENT[known_active] != unname(outcome_labels[alpha_arm[known_active]]))) ||
        ("REFERENCE" %in% supplied && "0" %in% names(outcome_labels) &&
         any(alpha$REFERENCE != unname(outcome_labels["0"])))) {
      stop("Treatment labels disagree between alpha and beta for ", analysis, "/", stratum, ".")
    }
  }
  matched <- match(alpha_key, beta_key)
  present <- !is.na(matched)
  diagnostics <- list(
    treatment_rows = nrow(alpha), outcome_rows = nrow(outcome),
    matched_rows = sum(present), treatment_only_rows = sum(!present),
    outcome_only_rows = sum(!beta_key %in% alpha_key),
    matched_outcome_rows = length(unique(matched[present]))
  )
  if (!any(present)) {
    stop("No matching ANALYTE_NAME x FU keys for ", analysis, "/", stratum, ".")
  }
  if (diagnostics$treatment_only_rows > 0L || diagnostics$outcome_only_rows > 0L) {
    warning("Partial mediation join for ", analysis, "/", stratum, ": ",
            diagnostics$matched_rows, " matched, ", diagnostics$treatment_only_rows,
            " treatment-only, and ", diagnostics$outcome_only_rows,
            " outcome-only rows. Treatment-only rows retain unavailable beta placeholders.")
  }
  mediation <- alpha
  mediation$BETA_EFFECT <- outcome$EFFECT_SIZE[matched]
  mediation$BETA_SE <- outcome$SE[matched]
  mediation$BETA_STATUS <- outcome$STATUS[matched]
  mediation$BETA_STATUS[!present] <- "not_estimable"
  mediation$STATUS <- ifelse(
    mediation$ALPHA_STATUS == "fit_error" | mediation$BETA_STATUS == "fit_error", "fit_error",
    ifelse(mediation$ALPHA_STATUS == "ok" & mediation$BETA_STATUS == "ok", "ok", "not_estimable")
  )
  if (treatment_filtered) {
    if (any(treatment$IN_FILTERED_SET[present] != outcome$IN_FILTERED_SET[matched[present]])) {
      stop("Filtered-probe membership disagrees for ", analysis, "/", stratum, ".")
    }
    mediation$IN_FILTERED_SET <- treatment$IN_FILTERED_SET
  }
  available <- mediation$STATUS == "ok"
  mediation$INDIRECT_EFFECT <- NA_real_
  mediation$INDIRECT_SE <- NA_real_
  mediation$INDIRECT_EFFECT[available] <- mediation$ALPHA_EFFECT[available] * mediation$BETA_EFFECT[available]
  mediation$INDIRECT_SE[available] <- sqrt(
    mediation$BETA_EFFECT[available]^2 * mediation$ALPHA_SE[available]^2 +
      mediation$ALPHA_EFFECT[available]^2 * mediation$BETA_SE[available]^2
  )
  mediation$Z_VALUE <- mediation$INDIRECT_EFFECT / mediation$INDIRECT_SE
  singular <- which(available & mediation$INDIRECT_SE == 0 & mediation$INDIRECT_EFFECT == 0)
  mediation$Z_VALUE[singular] <- 0
  mediation$P_VALUE <- 2 * pnorm(-abs(mediation$Z_VALUE))
  critical_value <- qnorm(0.975)
  mediation$CI_LOWER <- mediation$INDIRECT_EFFECT - critical_value * mediation$INDIRECT_SE
  mediation$CI_UPPER <- mediation$INDIRECT_EFFECT + critical_value * mediation$INDIRECT_SE
  mediation$BH_P_VALUE <- .adjust_within_fu(mediation$P_VALUE, mediation$FU, mediation$CONTRAST)
  if (treatment_filtered) {
    filtered_p <- mediation$P_VALUE
    filtered_p[!mediation$IN_FILTERED_SET] <- NA_real_
    mediation$BH_P_VALUE_FILTERED <- .adjust_within_fu(filtered_p, mediation$FU, mediation$CONTRAST)
  }
  if (outcome_type == "tte") mediation$INDIRECT_HAZARD_RATIO <- exp(mediation$INDIRECT_EFFECT)
  mediation <- mediation[order(mediation$ANALYTE_NAME, mediation$FU,
                               as.integer(sub("_vs_0$", "", mediation$CONTRAST))), , drop = FALSE]
  row.names(mediation) <- NULL
  list(mediation_effects = mediation, join_diagnostics = diagnostics)
}
